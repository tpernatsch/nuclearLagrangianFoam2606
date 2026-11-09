/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "DiffusiveMultiComponentPhaseModel.H"

#include "multiphaseInterSystem.H"
#include "fvcVolumeIntegrate.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvmLaplacian.H"
#include "fvcLaplacian.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcDDt.H"
#include "fvMatrix.H"
#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "MULES.H"
#include "CMULES.H"
#include "subCycle.H"
#include "IOdictionary.H"
#include "PstreamReduceOps.H"

namespace
{
    inline Foam::word baseSpecieNameForLookup(const Foam::word& s)
    {
        Foam::word m = Foam::IOobject::member(s);

        while (m.size() && (m[0] == 'Z' || m[0] == 'X'))
        {
            m = m.substr(1);
        }

        return m;
    }
}

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class BasePhaseModel, class phaseThermo>
Foam::DiffusiveMultiComponentPhaseModel<BasePhaseModel, phaseThermo>::
DiffusiveMultiComponentPhaseModel
(
    const multiphaseInterSystem& fluid,
    const word& phaseName
)
:
    BasePhaseModel(fluid, phaseName),
    species_(),
    inertIndex_(-1),
    thermoPtr_(nullptr),
    X_(),
    addDiffusion_(false),
    Dmol_("Dmol", dimArea/dimTime, 0.0),
    Sct_(1.0),
    Z_()
{
    thermoPtr_.reset
    (
        phaseThermo::New
        (
            fluid.mesh(),
            phaseName,
            basicThermo::phasePropertyName(basicThermo::dictName, phaseName)
        )
    );

    if (thermoPtr_->composition().species().empty())
    {
        FatalErrorInFunction
            << "The selected thermo is pure. Use a multicomponent thermo."
            << exit(FatalError);
    }

    species_ = thermoPtr_->composition().species();

    inertIndex_ = species_.find(thermoPtr_().template get<word>("inertSpecie"));

    if (inertIndex_ < 0)
    {
        FatalErrorInFunction
            << "Cannot find inertSpecie in species list for phase "
            << phaseName << exit(FatalError);
    }

    addDiffusion_ = thermoPtr_().template getOrDefault<bool>("addDiffusion", false);
    Sct_ = thermoPtr_().template getOrDefault<scalar>("Sct", 1.0);
    Dmol_.value() = thermoPtr_().template getOrDefault<scalar>("Dmol", 0.0);

    X_.setSize(species_.size());
    Z_.setSize(species_.size());

    boolList zRead(species_.size(), false);

    forAll(species_, i)
    {
        X_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("X" + species_[i], phaseName),
                    fluid.mesh().time().timeName(),
                    fluid.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                Y()[i]
            )
        );

        const word zName = IOobject::groupName("Z" + species_[i], phaseName);

        IOobject zHeader
        (
            zName,
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        );

        zRead[i] = zHeader.typeHeaderOk<volScalarField>(true);

        if (zRead[i])
        {
            Info<< "Reading restart field " << zName
                << " from time " << fluid.mesh().time().timeName() << nl;

            Z_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        zName,
                        fluid.mesh().time().timeName(),
                        fluid.mesh(),
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    fluid.mesh()
                )
            );
        }
        else
        {
            Info<< "Initializing field " << zName
                << " from alpha*Y and enabling AUTO_WRITE" << nl;

            Z_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        zName,
                        fluid.mesh().time().timeName(),
                        fluid.mesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    Y()[i]
                )
            );
        }

    }
        calculateVolumeFractions();

        forAll(species_, i)
        {
            if (!zRead[i])
            {
                Z_[i] = (*this)*Y()[i];
            }

            syncZBoundary(i);
            Z_[i].oldTime();
            fluid.mesh().setFluxRequired(Z_[i].name());
        }
}

template<class BasePhaseModel, class phaseThermo>
void Foam::DiffusiveMultiComponentPhaseModel<BasePhaseModel, phaseThermo>::
calculateVolumeFractions()
{
    volScalarField Xtotal(0.0*X_[0]);
    const volScalarField W(thermo().W());

    forAll(X_, i)
    {
        const dimensionedScalar Wi
        (
            "Wi",
            dimMass/dimMoles,
            thermo().composition().W(i)
        );

        if (i != inertIndex_)
        {
            X_[i] = W*Y()[i]/Wi;
            Xtotal += X_[i];
            X_[i].correctBoundaryConditions();

            const volScalarField::Boundary& YBf = Y()[i].boundaryField();

            forAll(YBf, patchi)
            {
                const fvPatchScalarField& YPf = YBf[patchi];

                if (YPf.fixesValue())
                {
                    scalarField& xbf = X_[i].boundaryFieldRef()[patchi];
                    const scalarField& ybf = Y()[i].boundaryField()[patchi];
                    const scalarField& Wbf = W.boundaryField()[patchi];

                    forAll(xbf, facei)
                    {
                        xbf[facei] = Wbf[facei]*ybf[facei]/Wi.value();
                    }
                }
            }
        }
    }

    X_[inertIndex_] = 1.0 - Xtotal;
    X_[inertIndex_].correctBoundaryConditions();
}


template<class BasePhaseModel, class phaseThermo>
void Foam::DiffusiveMultiComponentPhaseModel<BasePhaseModel, phaseThermo>::
calculateMassFractions()
{
    volScalarField W(X_[0]*thermo().composition().W(0));

    for (label i = 1; i < species_.size(); ++i)
    {
        W += X_[i]*thermo().composition().W(i);
    }

    forAll(Y(), i)
    {
        Y()[i] = X_[i]*thermo().composition().W(i)/W;

        Info<< Y()[i].name() << " mass fraction = "
            << " Min(Y)=" << min(Y()[i]).value()
            << " Max(Y)=" << max(Y()[i]).value()
            << nl;
    }
}


template<class BasePhaseModel, class phaseThermo>
const phaseThermo&
Foam::DiffusiveMultiComponentPhaseModel<BasePhaseModel, phaseThermo>::thermo() const
{
    return thermoPtr_();
}


template<class BasePhaseModel, class phaseThermo>
phaseThermo&
Foam::DiffusiveMultiComponentPhaseModel<BasePhaseModel, phaseThermo>::thermo()
{
    return thermoPtr_();
}


template<class BasePhaseModel, class phaseThermo>
void Foam::DiffusiveMultiComponentPhaseModel<BasePhaseModel, phaseThermo>::correct()
{
    thermo().correct();
}


template<class BasePhaseModel, class phaseThermo>
void Foam::DiffusiveMultiComponentPhaseModel<BasePhaseModel, phaseThermo>::syncZBoundary
(
    const label i
)
{
    const volScalarField& alpha1 = *this;
    const volScalarField& Yi = Y()[i];
    volScalarField& Zi = Z_[i];

    forAll(Zi.boundaryField(), patchi)
    {
        fvPatchScalarField& Zip = Zi.boundaryFieldRef()[patchi];
        const fvPatchScalarField& alphap = alpha1.boundaryField()[patchi];
        const fvPatchScalarField& Yip = Yi.boundaryField()[patchi];

        forAll(Zip, facei)
        {
            Zip[facei] = alphap[facei]*Yip[facei];
        }
    }

    Zi.correctBoundaryConditions();
}


template<class BasePhaseModel, class phaseThermo>
void Foam::DiffusiveMultiComponentPhaseModel<BasePhaseModel, phaseThermo>::syncZField
(
    const label i
)
{
    const volScalarField& alpha1 = *this;
    const volScalarField& Yi = Y()[i];
    volScalarField& Zi = Z_[i];

    scalarField& z = Zi.primitiveFieldRef();
    const scalarField& a = alpha1.primitiveField();
    const scalarField& y = Yi.primitiveField();

    forAll(z, celli)
    {
        z[celli] = a[celli]*y[celli];
    }

    syncZBoundary(i);
}


template<class BasePhaseModel, class phaseThermo>
void Foam::DiffusiveMultiComponentPhaseModel<BasePhaseModel, phaseThermo>::solveYi
(
    PtrList<volScalarField::Internal>& Su,
    PtrList<volScalarField::Internal>& Sp
)
{
    const volScalarField& alpha1 = *this;
    const fvMesh& mesh = alpha1.mesh();
    const scalarField& V = mesh.V();

    const dictionary& controls = mesh.solverDict(alpha1.name());

    const bool budgetYi = controls.getOrDefault<bool>("budgetYi", false);

    wordList budgetYiSpecies;
    if (budgetYi && controls.found("budgetYiSpecies"))
    {
        controls.lookup("budgetYiSpecies") >> budgetYiSpecies;
    }

    auto wantBudgetFor = [&](const word& name) -> bool
    {
        if (!budgetYi)
        {
            return false;
        }

        if (budgetYiSpecies.empty())
        {
            return true;
        }

        forAll(budgetYiSpecies, k)
        {
            if (budgetYiSpecies[k] == name)
            {
                return true;
            }
        }

        return false;
    };

    auto alphaCutoffForPhase = [](const fvMesh& mesh, const word& phName) -> scalar
    {
        const scalar defCut = 1e-6;
        scalar cut = defCut;

        if (mesh.foundObject<IOdictionary>("phaseProperties"))
        {
            const IOdictionary& pp = mesh.lookupObject<IOdictionary>("phaseProperties");

            if (pp.found(phName))
            {
                const dictionary& pd = pp.subDict(phName);
                cut = pd.getOrDefault<scalar>
                (
                    "residualAlpha",
                    pd.getOrDefault<scalar>("alphaCutoff", defCut)
                );
            }
            else
            {
                cut = pp.getOrDefault<scalar>
                (
                    "residualAlpha",
                    pp.getOrDefault<scalar>("alphaCutoff", defCut)
                );
            }
        }
        else
        {
            IOdictionary pp
            (
                IOobject
                (
                    "phaseProperties",
                    mesh.time().constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            if (pp.found(phName))
            {
                const dictionary& pd = pp.subDict(phName);
                cut = pd.getOrDefault<scalar>
                (
                    "residualAlpha",
                    pd.getOrDefault<scalar>("alphaCutoff", defCut)
                );
            }
            else
            {
                cut = pp.getOrDefault<scalar>
                (
                    "residualAlpha",
                    pp.getOrDefault<scalar>("alphaCutoff", defCut)
                );
            }
        }

        return cut > 0 ? max(cut, SMALL) : scalar(0);
    };

    const scalar alphaCut = alphaCutoffForPhase(mesh, this->name());
    const scalar alphaSolveCut = max(alphaCut, SMALL);

    const bool cleanupGhostYi = controls.getOrDefault<bool>("cleanupGhostYi", true);

    const scalar alphaWipe =
        controls.getOrDefault<scalar>("alphaCleanupCut", alphaSolveCut) > 0
      ? max(controls.getOrDefault<scalar>("alphaCleanupCut", alphaSolveCut), SMALL)
      : alphaSolveCut;

    const scalar alphaSpeciesCut =
        controls.getOrDefault<scalar>("alphaSpeciesCut", alphaWipe) > 0
      ? max(controls.getOrDefault<scalar>("alphaSpeciesCut", alphaWipe), SMALL)
      : alphaWipe;

    const bool supportProjectBeforeSolve =
        controls.getOrDefault<bool>("supportProjectBeforeSolve", true);

    const scalar ZboundMassTol =
        controls.getOrDefault<scalar>("ZboundMassTol", 1e-30);

    const label nSub =
        max(label(1), controls.getOrDefault<label>("nYiSubCycles", 1));

    const scalar globalDt = mesh.time().deltaTValue();
    const scalar yiDt = globalDt/scalar(nSub);

    const bool printYiSetup = (nSub > 1 || addDiffusion_);

    if (printYiSetup)
    {
        static HashTable<label> printedYiSetup;

        const word setupKey = this->name();

        if (!printedYiSetup.found(setupKey))
        {
            if (Pstream::master())
            {
                Info<< "YI Zimplicit phase=" << this->name()
                    << " nYiSubCycles=" << nSub
                    << " globalDt=" << globalDt
                    << " yiDt=" << yiDt
                    << " addDiffusion=" << addDiffusion_
                    << " Dmol=" << Dmol_.value()
                    << " Sct=" << Sct_
                    << nl;
            }

            printedYiSetup.set(setupKey, 1);
        }
    }

    surfaceScalarField phiPhase
    (
        IOobject
        (
            "phiPhase_" + this->name(),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->alphaPhi()
    );
    phiPhase.setOriented(false);

    auto reconstructYFromZ = [&](const label i)
    {
        const volScalarField& Zi = Z_[i];
        volScalarField& Yi = Y()[i];

        forAll(Yi, celli)
        {
            if (alpha1[celli] <= alphaSpeciesCut)
            {
                Yi[celli] = 0.0;
            }
            else
            {
                scalar y = Zi[celli]/max(alpha1[celli], SMALL);

                if (y != y)
                {
                    y = 0.0;
                }

                Yi[celli] = min(max(y, scalar(0)), scalar(1));
            }
        }

        forAll(Yi.boundaryField(), patchi)
        {
            fvPatchScalarField& Yip = Yi.boundaryFieldRef()[patchi];
            const fvPatchScalarField& alphap = alpha1.boundaryField()[patchi];
            const fvPatchScalarField& Zip = Zi.boundaryField()[patchi];

            forAll(Yip, facei)
            {
                if (alphap[facei] <= alphaSpeciesCut)
                {
                    Yip[facei] = 0.0;
                }
                else
                {
                    scalar y = Zip[facei]/max(alphap[facei], SMALL);

                    if (y != y)
                    {
                        y = 0.0;
                    }

                    Yip[facei] = min(max(y, scalar(0)), scalar(1));
                }
            }
        }

        Yi.correctBoundaryConditions();
    };

    auto conservativeBoundZ =
    [&]
    (
        volScalarField& Zi,
        const volScalarField& rho1
    ) -> scalar
    {
        scalarField& z = Zi.primitiveFieldRef();

        scalar dMClampLocal = 0.0;
        scalar capMassLocal = 0.0;
        scalar stockMassLocal = 0.0;
        scalar movedMassLocal = 0.0;

        forAll(z, celli)
        {
            const scalar a = max(alpha1[celli], scalar(0));
            const scalar zOld = z[celli];
            const scalar zNew = min(max(zOld, scalar(0)), a);

            if (zNew != zOld)
            {
                dMClampLocal += rho1[celli]*V[celli]*(zNew - zOld);
                z[celli] = zNew;
            }
        }

        scalar dMClamp = dMClampLocal;
        reduce(dMClamp, sumOp<scalar>());

        if (mag(dMClamp) > ZboundMassTol)
        {
            if (dMClamp < 0.0)
            {
                forAll(z, celli)
                {
                    if (alpha1[celli] > alphaSpeciesCut)
                    {
                        capMassLocal +=
                            rho1[celli]*V[celli]
                           *max(alpha1[celli] - z[celli], scalar(0));
                    }
                }

                scalar capMass = capMassLocal;
                reduce(capMass, sumOp<scalar>());

                if (capMass > ZboundMassTol)
                {
                    const scalar frac = min((-dMClamp)/capMass, scalar(1));

                    forAll(z, celli)
                    {
                        if (alpha1[celli] > alphaSpeciesCut)
                        {
                            const scalar dz =
                                frac*max(alpha1[celli] - z[celli], scalar(0));

                            z[celli] += dz;
                            movedMassLocal += rho1[celli]*V[celli]*dz;
                        }
                    }
                }
            }
            else
            {
                forAll(z, celli)
                {
                    if (alpha1[celli] > alphaSpeciesCut)
                    {
                        stockMassLocal += rho1[celli]*V[celli]*max(z[celli], scalar(0));
                    }
                }

                scalar stockMass = stockMassLocal;
                reduce(stockMass, sumOp<scalar>());

                if (stockMass > ZboundMassTol)
                {
                    const scalar frac = min(dMClamp/stockMass, scalar(1));

                    forAll(z, celli)
                    {
                        if (alpha1[celli] > alphaSpeciesCut)
                        {
                            const scalar dz = frac*z[celli];
                            z[celli] -= dz;
                            movedMassLocal += rho1[celli]*V[celli]*dz;
                        }
                    }
                }
            }
        }

        scalar movedMass = movedMassLocal;
        reduce(movedMass, sumOp<scalar>());

        return dMClamp < 0.0 ? dMClamp + movedMass : dMClamp - movedMass;
    };

    auto makeGammaDiff = [&]() -> autoPtr<surfaceScalarField>
    {
        autoPtr<surfaceScalarField> GammaDiffPtr;

        if (!addDiffusion_ || Dmol_.value() <= 0.0)
        {
            return GammaDiffPtr;
        }

        scalar alphaDiffMask =
            controls.getOrDefault<scalar>("alphaDiffMask", alphaSpeciesCut);

        alphaDiffMask = max(alphaDiffMask, scalar(0));

        surfaceScalarField alphaDiffF
        (
            IOobject
            (
                "alphaDiffF_" + this->name(),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0)
        );

        {
            scalarField& aF = alphaDiffF.primitiveFieldRef();
            const labelUList& own = mesh.owner();
            const labelUList& nei = mesh.neighbour();

            forAll(aF, facei)
            {
                const scalar aOwn = max(alpha1[own[facei]], scalar(0));
                const scalar aNei = max(alpha1[nei[facei]], scalar(0));

                aF[facei] =
                    (aOwn > alphaDiffMask && aNei > alphaDiffMask)
                  ? min(aOwn, aNei)
                  : scalar(0);
            }
        }

        forAll(alphaDiffF.boundaryField(), patchi)
        {
            fvsPatchScalarField& aPf = alphaDiffF.boundaryFieldRef()[patchi];
            const fvPatchScalarField& aOwnPf = alpha1.boundaryField()[patchi];

            if (aPf.coupled())
            {
                tmp<scalarField> taNbr = aOwnPf.patchNeighbourField();
                const scalarField& aNbr = taNbr();

                forAll(aPf, facei)
                {
                    const scalar aOwn = max(aOwnPf[facei], scalar(0));
                    const scalar aNei = max(aNbr[facei], scalar(0));

                    aPf[facei] =
                        (aOwn > alphaDiffMask && aNei > alphaDiffMask)
                      ? min(aOwn, aNei)
                      : scalar(0);
                }
            }
            else
            {
                forAll(aPf, facei)
                {
                    const scalar aOwn = max(aOwnPf[facei], scalar(0));
                    aPf[facei] = aOwn > alphaDiffMask ? aOwn : scalar(0);
                }
            }
        }

        const surfaceScalarField nutF(fvc::interpolate(this->fluid().turbulence()->nut()));

        GammaDiffPtr.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "GammaDiffF_" + this->name(),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                alphaDiffF*(Dmol_ + nutF/max(Sct_, SMALL))
            )
        );

        return GammaDiffPtr;
    };

    for (label subI = 0; subI < nSub; ++subI)
    {
        const scalar dt = yiDt;
        const dimensionedScalar rDeltaT("rDeltaT", dimless/dimTime, 1.0/dt);

        const tmp<volScalarField> trho = this->rho();
        const volScalarField& rho1 = trho();

        autoPtr<surfaceScalarField> GammaDiffPtr = makeGammaDiff();

        volScalarField Ysum(0.0*Y()[0]);

        forAll(Y(), i)
        {
            if (i == inertIndex_)
            {
                continue;
            }

            volScalarField& Yi = Y()[i];
            volScalarField& Zi = Z_[i];

            const bool doBudget = wantBudgetFor(Yi.name());
            const scalar M0 =
                doBudget ? fvc::domainIntegrate(rho1*Zi).value() : scalar(0);

            syncZBoundary(i);

            if (supportProjectBeforeSolve)
            {
                conservativeBoundZ(Zi, rho1);
                syncZBoundary(i);
            }

            reconstructYFromZ(i);

            const scalarField zStart(Zi.primitiveField());

            const word divSchemeYi("div(phi," + Yi.name() + ')');
            const dictionary& YiSolverDict = mesh.solverDict(Yi.name());

            const dimensionedScalar alphaFloor
            (
                "alphaSolveCut",
                alpha1.dimensions(),
                alphaSolveCut
            );

            const tmp<volScalarField> tAlphaSolve(max(alpha1, alphaFloor));
            const volScalarField& alphaSolve = tAlphaSolve();

            surfaceScalarField alphaFaceConv
            (
                IOobject
                (
                    "alphaFaceConv_" + Yi.name(),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                max
                (
                    fvc::interpolate(alphaSolve),
                    dimensionedScalar("alphaFaceFloor", dimless, alphaSolveCut)
                )
            );

            const tmp<surfaceScalarField> tPhiAlphaByScheme =
                fvc::flux(phiPhase, alphaSolve, divSchemeYi);

            const surfaceScalarField& phiAlphaByScheme = tPhiAlphaByScheme();

            {
                scalarField& aF = alphaFaceConv.primitiveFieldRef();
                const scalarField& phiF = phiPhase.primitiveField();
                const scalarField& paF = phiAlphaByScheme.primitiveField();

                forAll(aF, facei)
                {
                    if (mag(phiF[facei]) > VSMALL)
                    {
                        aF[facei] = paF[facei]/phiF[facei];
                    }

                    aF[facei] = min(max(aF[facei], alphaSolveCut), scalar(1));
                }
            }

            forAll(alphaFaceConv.boundaryField(), patchi)
            {
                fvsPatchScalarField& aFp = alphaFaceConv.boundaryFieldRef()[patchi];
                const fvsPatchScalarField& phiP = phiPhase.boundaryField()[patchi];
                const fvsPatchScalarField& paP = phiAlphaByScheme.boundaryField()[patchi];

                forAll(aFp, facei)
                {
                    if (mag(phiP[facei]) > VSMALL)
                    {
                        aFp[facei] = paP[facei]/phiP[facei];
                    }

                    aFp[facei] = min(max(aFp[facei], alphaSolveCut), scalar(1));
                }
            }

            surfaceScalarField phiOverAlpha
            (
                IOobject
                (
                    "phiOverAlpha_" + Yi.name(),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                phiPhase/alphaFaceConv
            );

            volScalarField::Internal SuA
            (
                IOobject
                (
                    "SuA_" + Yi.name(),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensionedScalar("zero", Su[i].dimensions(), 0.0)
            );

            SuA.field() = alpha1.primitiveField()*Su[i].field();

            volScalarField::Internal DiffA
            (
                IOobject
                (
                    "DiffA_" + Yi.name(),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensionedScalar("zero", dimless/dimTime, 0.0)
            );

            if (GammaDiffPtr.valid())
            {
                surfaceScalarField phiDiffY
                (
                    IOobject
                    (
                        "phiDiffY_" + Yi.name(),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    (*GammaDiffPtr)*fvc::snGrad(Yi)*mesh.magSf()
                );

                const tmp<volScalarField> tDivPhiDiffY = fvc::div(phiDiffY);
                const volScalarField& divPhiDiffY = tDivPhiDiffY();

                DiffA.field() = divPhiDiffY.primitiveField();

                if (controls.getOrDefault<bool>("conserveFaceFluxY", true))
                {
                    const tmp<surfaceScalarField> trhof = fvc::interpolate(rho1);
                    const surfaceScalarField& rhof = trhof();

                    scalar physicalPatchMassRateLocal = 0.0;

                    forAll(mesh.boundary(), patchi)
                    {
                        const fvPatch& p = mesh.boundary()[patchi];

                        if (!p.coupled())
                        {
                            physicalPatchMassRateLocal +=
                                sum
                                (
                                    rhof.boundaryField()[patchi]
                                   *phiDiffY.boundaryField()[patchi]
                                );
                        }
                    }

                    scalar physicalPatchMassRate = physicalPatchMassRateLocal;
                    reduce(physicalPatchMassRate, sumOp<scalar>());

                    scalar rawDiffMassRateLocal = 0.0;
                    scalar activeMassLocal = 0.0;

                    forAll(DiffA, celli)
                    {
                        rawDiffMassRateLocal += rho1[celli]*DiffA[celli]*V[celli];

                        if (alpha1[celli] > alphaSpeciesCut)
                        {
                            activeMassLocal += rho1[celli]*V[celli];
                        }
                    }

                    scalar rawDiffMassRate = rawDiffMassRateLocal;
                    scalar activeMass = activeMassLocal;

                    reduce(rawDiffMassRate, sumOp<scalar>());
                    reduce(activeMass, sumOp<scalar>());

                    if (activeMass > VSMALL)
                    {
                        const scalar correctionA =
                            (rawDiffMassRate - physicalPatchMassRate)/activeMass;

                        forAll(DiffA, celli)
                        {
                            if (alpha1[celli] > alphaSpeciesCut)
                            {
                                DiffA[celli] -= correctionA;
                            }
                        }
                    }
                }
            }

            volScalarField::Internal ZDdtStartA
            (
                IOobject
                (
                    "ZDdtStartA_" + Yi.name(),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensionedScalar("zero", dimless/dimTime, 0.0)
            );

            forAll(ZDdtStartA, celli)
            {
                ZDdtStartA[celli] = zStart[celli]/dt;
            }

            fvScalarMatrix ZiEqn
            (
                fvm::Sp(rDeltaT, Zi)
              + fvm::div(phiOverAlpha, Zi, divSchemeYi)
             ==
                ZDdtStartA
              + SuA
              + DiffA
              + fvm::Sp(Sp[i], Zi)
            );

            ZiEqn.solve(YiSolverDict);

            conservativeBoundZ(Zi, rho1);
            syncZBoundary(i);
            reconstructYFromZ(i);

            if (cleanupGhostYi)
            {
                forAll(Yi, celli)
                {
                    if (alpha1[celli] <= alphaWipe)
                    {
                        Yi[celli] = 0.0;
                    }
                }
            }

            Yi.correctBoundaryConditions();
            Ysum += Yi;

            if (doBudget)
            {
                const scalar M1 = fvc::domainIntegrate(rho1*Zi).value();

                scalar MsourceLocal = 0.0;
                forAll(SuA, celli)
                {
                    MsourceLocal += rho1[celli]*SuA[celli]*V[celli]*dt;
                }

                scalar Msource = MsourceLocal;
                reduce(Msource, sumOp<scalar>());

                Info<< "YI BUDGET " << Yi.name()
                    << " phase=" << this->name()
                    << " subStep=" << subI + 1 << '/' << nSub
                    << " M0=" << M0
                    << " M1=" << M1
                    << " dM=" << M1 - M0
                    << " Msource=" << Msource
                    << " nonSource=" << (M1 - M0) - Msource
                    << nl;
            }
        }

        forAll(Ysum, celli)
        {
            if (Ysum[celli] > 1.0 + 1e-12)
            {
                const scalar inv = 1.0/Ysum[celli];

                forAll(Y(), j)
                {
                    if (j != inertIndex_)
                    {
                        Y()[j][celli] *= inv;
                    }
                }

                Ysum[celli] = 1.0;
            }
        }

        volScalarField& Yinert = Y()[inertIndex_];
        Yinert = scalar(1) - Ysum;
        Yinert = max(Yinert, dimensionedScalar("zero", Yinert.dimensions(), 0.0));

        if (cleanupGhostYi)
        {
            forAll(Yinert, celli)
            {
                if (alpha1[celli] <= alphaWipe)
                {
                    Yinert[celli] = 1.0;
                }
            }
        }

        Yinert.correctBoundaryConditions();

        volScalarField& Zinert = Z_[inertIndex_];

        forAll(alpha1, celli)
        {
            const scalar a = max(alpha1[celli], scalar(0));
            scalar sumZNonInert = 0.0;

            forAll(Y(), j)
            {
                if (j != inertIndex_)
                {
                    sumZNonInert += max(Z_[j][celli], scalar(0));
                }
            }

            Zinert[celli] = max(min(a - sumZNonInert, a), scalar(0));
        }

        forAll(Y(), i)
        {
            syncZBoundary(i);
        }
    }

    forAll(Y(), i)
    {
        Info<< Y()[i].name() << " mass fraction = "
            << " Min(Y)=" << min(Y()[i]).value()
            << " Max(Y)=" << max(Y()[i]).value()
            << nl;
    }
}


template<class BasePhaseModel, class phaseThermo>
const Foam::PtrList<Foam::volScalarField>&
Foam::DiffusiveMultiComponentPhaseModel<BasePhaseModel, phaseThermo>::Y() const
{
    return thermoPtr_->composition().Y();
}


template<class BasePhaseModel, class phaseThermo>
Foam::PtrList<Foam::volScalarField>&
Foam::DiffusiveMultiComponentPhaseModel<BasePhaseModel, phaseThermo>::Y()
{
    return thermoPtr_->composition().Y();
}


template<class BasePhaseModel, class phaseThermo>
Foam::label
Foam::DiffusiveMultiComponentPhaseModel<BasePhaseModel, phaseThermo>::inertIndex() const
{
    return inertIndex_;
}


