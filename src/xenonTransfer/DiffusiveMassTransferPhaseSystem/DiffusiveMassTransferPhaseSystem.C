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

#include "DiffusiveMassTransferPhaseSystem.H"
#include "HashPtrTable.H"
#include "fvcDiv.H"
#include "fvmSup.H"
#include "fvcGrad.H"
#include "fvMatrix.H"
#include "volFields.H"
#include "fundamentalConstants.H"
#include "PstreamReduceOps.H"

using namespace Foam::multiphaseInter;

namespace
{
    inline Foam::word baseSpecieName(const Foam::word& s)
    {
        Foam::word m = Foam::IOobject::member(s);

        while (m.size() && (m[0] == 'Z' || m[0] == 'X'))
        {
            m = m.substr(1);
        }

        return m;
    }


    inline Foam::volScalarField* volScalarPtr
    (
        const Foam::fvMesh& mesh,
        const Foam::word& name
    )
    {
        if (mesh.foundObject<Foam::volScalarField>(name))
        {
            return const_cast<Foam::volScalarField*>
            (
                &mesh.lookupObject<Foam::volScalarField>(name)
            );
        }

        return nullptr;
    }


    inline bool syncSpeciesStorage
    (
        const Foam::phaseModel& phase,
        const Foam::word& specie,
        const Foam::scalar alphaCut,
        const Foam::scalar ZboundMassTol
    )
    {
        const Foam::volScalarField& alpha = phase;
        const Foam::fvMesh& mesh = alpha.mesh();
        const Foam::scalarField& V = mesh.V();

        const Foam::word phaseName = phase.name();

        const Foam::word Yname =
            Foam::IOobject::groupName(specie, phaseName);

        const Foam::word Zname =
            Foam::IOobject::groupName("Z" + specie, phaseName);

        Foam::volScalarField* Yptr = volScalarPtr(mesh, Yname);
        Foam::volScalarField* Zptr = volScalarPtr(mesh, Zname);

        if (!Yptr || !Zptr)
        {
            return false;
        }

        Foam::volScalarField& Y = *Yptr;
        Foam::volScalarField& Z = *Zptr;

        const Foam::tmp<Foam::volScalarField> trho = phase.rho();
        const Foam::volScalarField& rho = trho();

        Foam::scalarField& z = Z.primitiveFieldRef();

        Foam::scalar dMClampLocal = 0.0;
        Foam::scalar capMassLocal = 0.0;
        Foam::scalar stockMassLocal = 0.0;
        Foam::scalar movedMassLocal = 0.0;

        Foam::label nHighLocal = 0;
        Foam::label nLowLocal = 0;

        forAll(z, celli)
        {
            const Foam::scalar a = Foam::max(alpha[celli], Foam::scalar(0));
            const Foam::scalar zOld = z[celli];
            const Foam::scalar zNew = Foam::min(Foam::max(zOld, Foam::scalar(0)), a);

            if (zOld > a + 1e-12)
            {
                ++nHighLocal;
            }

            if (zOld < -1e-12)
            {
                ++nLowLocal;
            }

            if (zNew != zOld)
            {
                dMClampLocal += rho[celli]*V[celli]*(zNew - zOld);
                z[celli] = zNew;
            }
        }

        Foam::scalar dMClamp = dMClampLocal;
        Foam::label nHigh = nHighLocal;
        Foam::label nLow = nLowLocal;

        Foam::reduce(dMClamp, Foam::sumOp<Foam::scalar>());
        Foam::reduce(nHigh, Foam::sumOp<Foam::label>());
        Foam::reduce(nLow, Foam::sumOp<Foam::label>());

        // Redistribute any clamped mass (see syncAndSnapshotSpecies)
        if (Foam::mag(dMClamp) > ZboundMassTol)
        {
            if (dMClamp < 0.0)
            {
                forAll(z, celli)
                {
                    if (alpha[celli] > alphaCut)
                    {
                        const Foam::scalar cap =
                            Foam::max(alpha[celli] - z[celli], Foam::scalar(0));

                        capMassLocal += rho[celli]*V[celli]*cap;
                    }
                }

                Foam::scalar capMass = capMassLocal;
                Foam::reduce(capMass, Foam::sumOp<Foam::scalar>());

                if (capMass > ZboundMassTol)
                {
                    const Foam::scalar frac =
                        Foam::min((-dMClamp)/capMass, Foam::scalar(1));

                    forAll(z, celli)
                    {
                        if (alpha[celli] > alphaCut)
                        {
                            const Foam::scalar cap =
                                Foam::max(alpha[celli] - z[celli], Foam::scalar(0));

                            const Foam::scalar dz = frac*cap;

                            z[celli] += dz;
                            movedMassLocal += rho[celli]*V[celli]*dz;
                        }
                    }
                }
            }
            else
            {
                forAll(z, celli)
                {
                    if (alpha[celli] > alphaCut)
                    {
                        stockMassLocal +=
                            rho[celli]*V[celli]*Foam::max(z[celli], Foam::scalar(0));
                    }
                }

                Foam::scalar stockMass = stockMassLocal;
                Foam::reduce(stockMass, Foam::sumOp<Foam::scalar>());

                if (stockMass > ZboundMassTol)
                {
                    const Foam::scalar frac =
                        Foam::min(dMClamp/stockMass, Foam::scalar(1));

                    forAll(z, celli)
                    {
                        if (alpha[celli] > alphaCut)
                        {
                            const Foam::scalar dz =
                                frac*Foam::max(z[celli], Foam::scalar(0));

                            z[celli] -= dz;
                            movedMassLocal += rho[celli]*V[celli]*dz;
                        }
                    }
                }
            }
        }

        Foam::scalar movedMass = movedMassLocal;
        Foam::reduce(movedMass, Foam::sumOp<Foam::scalar>());

        Foam::scalarField& y = Y.primitiveFieldRef();

        forAll(y, celli)
        {
            if (alpha[celli] <= alphaCut)
            {
                y[celli] = 0.0;
            }
            else
            {
                Foam::scalar yi = z[celli]/Foam::max(alpha[celli], Foam::SMALL);

                if (yi != yi)
                {
                    yi = 0.0;
                }

                y[celli] = Foam::min(Foam::max(yi, Foam::scalar(0)), Foam::scalar(1));
            }
        }

        Y.correctBoundaryConditions();

        forAll(Z.boundaryField(), patchi)
        {
            Foam::fvPatchScalarField& Zp = Z.boundaryFieldRef()[patchi];

            const Foam::fvPatchScalarField& ap =
                alpha.boundaryField()[patchi];

            const Foam::fvPatchScalarField& Yp =
                Y.boundaryField()[patchi];

            forAll(Zp, facei)
            {
                Zp[facei] = ap[facei]*Yp[facei];
            }
        }

        Z.correctBoundaryConditions();

        return true;
    }
    
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::DiffusiveMassTransferPhaseSystem<BasePhaseSystem>::DiffusiveMassTransferPhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh)
{
    this->generatePairsAndSubModels("massTransferModel", massTransferModels_);

    forAllConstIters(massTransferModels_, iterModel)
    {
        if (!dmdt_.found(iterModel()->pair()))
        {
            dmdt_.set
            (
                iterModel()->pair(),
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName("dmdt",iterModel()->pair().name()),
                        this->mesh().time().timeName(),
                        this->mesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    this->mesh(),
                    dimensionedScalar(dimDensity/dimTime, Zero)
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::DiffusiveMassTransferPhaseSystem<BasePhaseSystem>::calculateL
(
    const volScalarField& dmdtNetki,
    const phasePairKey& keyik,
    const phasePairKey& keyki,
    const volScalarField& T
) const
{
    auto tL = tmp<volScalarField>::New
    (
        IOobject
        (
            "tL",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimEnergy/dimMass, Zero)
    );
    auto& L = tL.ref();

    if (massTransferModels_.found(keyik))
    {
        const autoPtr<interfaceCompositionModel>& interfacePtr =
            massTransferModels_[keyik];

        word speciesName = interfacePtr->transferSpecie();

        const word species(speciesName.substr(0, speciesName.find('.')));

        L -= neg(dmdtNetki)*interfacePtr->L(species, T);
    }

    if (massTransferModels_.found(keyki))
    {
        const autoPtr<interfaceCompositionModel>& interfacePtr =
            massTransferModels_[keyki];

        word speciesName = interfacePtr->transferSpecie();

        const word species(speciesName.substr(0, speciesName.find('.')));

        L -= pos(dmdtNetki)*interfacePtr->L(species, T);
    }

    return tL;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::DiffusiveMassTransferPhaseSystem<BasePhaseSystem>::dmdt
(
    const phasePairKey& key
) const
{
    auto tdmdt = tmp<volScalarField>::New
    (
        IOobject
        (
            "dmdt",
            this->mesh().time().timeName(),
            this->mesh()
        ),
        this->mesh(),
        dimensionedScalar(dimDensity/dimTime, Zero)
    );
    auto& dmdt = tdmdt.ref();

    if (dmdt_.found(key))
    {
        dmdt = *dmdt_[key];
    }

    return tdmdt;
}


template<class BasePhaseSystem>
Foam::tmp<Foam::fvScalarMatrix>
Foam::DiffusiveMassTransferPhaseSystem<BasePhaseSystem>::heatTransfer
(
    const volScalarField& T
)
{
    auto teqn = tmp<fvScalarMatrix>::New(T, dimEnergy/dimTime);
    auto& eqn = teqn.ref();

    forAllConstIters(this->phaseModels_, iteri)
    {
        const phaseModel& phasei = iteri()();

        auto iterk = iteri;

        for (++iterk; iterk != this->phaseModels_.end(); ++iterk)
        {
            if (iteri()().name() != iterk()().name())
            {
                const phaseModel& phasek = iterk()();

                const phasePairKey keyik(phasei.name(), phasek.name(), true);

                const phasePairKey keyki(phasek.name(), phasei.name(), true);

                auto tdmdtNetki = tmp<volScalarField>::New
                (
                    IOobject
                    (
                        "tdmdtYki",
                        this->mesh().time().timeName(),
                        this->mesh()
                    ),
                    this->mesh(),
                    dimensionedScalar(dimDensity/dimTime, Zero)
                );
                auto& dmdtNetki = tdmdtNetki.ref();

                auto tSp = tmp<volScalarField>::New
                (
                    IOobject
                    (
                        "Sp",
                        this->mesh().time().timeName(),
                        this->mesh()
                    ),
                    this->mesh(),
                    dimensionedScalar(dimDensity/dimTime/dimTemperature, Zero)
                );
                auto& Sp = tSp.ref();

                auto tSu = tmp<volScalarField>::New
                (
                    IOobject
                    (
                        "Su",
                        this->mesh().time().timeName(),
                        this->mesh()
                    ),
                    this->mesh(),
                    dimensionedScalar(dimDensity/dimTime, Zero)
                );
                auto& Su = tSu.ref();


                if (massTransferModels_.found(keyik))
                {
                    autoPtr<interfaceCompositionModel>& interfacePtr =
                        massTransferModels_[keyik];

                    dmdtNetki -= *dmdt_[keyik];

                    tmp<volScalarField> KSp =
                        interfacePtr->KSp(interfaceCompositionModel::T, T);

                    if (KSp.valid())
                    {
                        Sp += KSp.ref();
                    }

                    tmp<volScalarField> KSu =
                        interfacePtr->KSu(interfaceCompositionModel::T, T);

                    if (KSu.valid())
                    {
                        Su += KSu.ref();
                    }

                    if (!KSp.valid() && !KSu.valid())
                    {
                        Su += *dmdt_[keyik];
                    }
                }

                if (massTransferModels_.found(keyki))
                {
                    autoPtr<interfaceCompositionModel>& interfacePtr =
                        massTransferModels_[keyki];

                    dmdtNetki += *dmdt_[keyki];


                    tmp<volScalarField> KSp =
                        interfacePtr->KSp(interfaceCompositionModel::T, T);

                    if (KSp.valid())
                    {
                        Sp += KSp.ref();
                    }

                    tmp<volScalarField> KSu =
                        interfacePtr->KSu(interfaceCompositionModel::T, T);

                    if (KSu.valid())
                    {
                        Su += KSu.ref();
                    }

                    if (!KSp.valid() && !KSu.valid())
                    {
                        Su += *dmdt_[keyki];
                    }
                }

                tmp<volScalarField> L = calculateL(dmdtNetki, keyik, keyki, T);

                eqn -= fvm::Sp(Sp*L.ref(), T) + Su*L.ref();
            }
        }
    }
    return teqn;
}


template<class BasePhaseSystem>
Foam::tmp<Foam::fvScalarMatrix>
Foam::DiffusiveMassTransferPhaseSystem<BasePhaseSystem>::volTransfer
(
    const volScalarField& p
)
{
    auto teqn = tmp<fvScalarMatrix>::New(p, dimVolume/dimTime);
    auto& eqn = teqn.ref();

    auto tSp = tmp<volScalarField>::New
    (
        IOobject
        (
            "Sp",
            this->mesh().time().timeName(),
            this->mesh()
        ),
        this->mesh(),
        dimensionedScalar(dimless/dimTime/dimPressure, Zero)
    );
    auto& Sp = tSp.ref();

    auto tSu = tmp<volScalarField>::New
    (
        IOobject
        (
            "Su",
            this->mesh().time().timeName(),
            this->mesh()
        ),
        this->mesh(),
        dimensionedScalar(dimless/dimTime, Zero)
    );
    auto& Su = tSu.ref();

    forAllConstIters(this->totalPhasePairs(), iter)
    {
        const phasePair& pair = iter()();

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();

        const phasePairKey key12
        (
            phase1.name(),
            phase2.name(),
            true
        );

        if (massTransferModels_.found(key12))
        {
            autoPtr<interfaceCompositionModel>& interfacePtr =
                massTransferModels_[key12];

            tmp<volScalarField> KSp =
                interfacePtr->KSp(interfaceCompositionModel::P, p);

            if (KSp.valid())
            {
                Sp -=
                    KSp.ref()
                   *(
                        - this->coeffs(phase1.name())
                        + this->coeffs(phase2.name())
                    );
            }

            tmp<volScalarField> KSu =
                interfacePtr->KSu(interfaceCompositionModel::P, p);

            if (KSu.valid())
            {
                Su -=
                    KSu.ref()
                   *(
                        - this->coeffs(phase1.name())
                        + this->coeffs(phase2.name())
                    );
            }

            if (!KSp.valid() && !KSu.valid())
            {
                Su -=
                    *dmdt_[key12]
                    *(
                        - this->coeffs(phase1.name())
                        + this->coeffs(phase2.name())
                    );
            }
        }

        const phasePairKey key21
        (
            phase2.name(),
            phase1.name(),
            true
        );

        if (massTransferModels_.found(key21))
        {
            autoPtr<interfaceCompositionModel>& interfacePtr =
                massTransferModels_[key21];

            tmp<volScalarField> KSp =
                interfacePtr->KSp(interfaceCompositionModel::P, p);

            if (KSp.valid())
            {
                Sp +=
                    KSp.ref()
                   *(
                        - this->coeffs(phase1.name())
                        + this->coeffs(phase2.name())
                    );
            }

            tmp<volScalarField> KSu =
                interfacePtr->KSu(interfaceCompositionModel::P, p);

            if (KSu.valid())
            {
                Su +=
                    KSu.ref()
                   *(
                        - this->coeffs(phase1.name())
                        + this->coeffs(phase2.name())
                    );
            }

            if (!KSp.valid() && !KSu.valid())
            {
                Su +=
                    *dmdt_[key21]
                    *(
                        - this->coeffs(phase1.name())
                        + this->coeffs(phase2.name())
                    );
            }
        }

    }

    eqn += fvm::Sp(Sp, p) + Su;
    return teqn;
}


template<class BasePhaseSystem>
void Foam::DiffusiveMassTransferPhaseSystem<BasePhaseSystem>::correctMassSources
(
    const volScalarField& T
)
{
    forAllConstIters(massTransferModels_, iter)
    {
        autoPtr<interfaceCompositionModel>& interfacePtr =
            massTransferModels_[iter.key()];

        const phasePair& pair = interfacePtr->pair();

        const word specie = baseSpecieName(interfacePtr->transferSpecie());

        const word alphaControlsName =
            IOobject::groupName("alpha", pair.to().name());

        const dictionary& alphaControls =
            this->mesh().solverDict(alphaControlsName);

        const scalar alphaSpeciesCut =
            max
            (
                alphaControls.getOrDefault<scalar>("alphaSpeciesCut", 1e-12),
                SMALL
            );

        const scalar ZboundMassTol =
            alphaControls.getOrDefault<scalar>("ZboundMassTol", 1e-30);

        syncSpeciesStorage(pair.from(), specie, alphaSpeciesCut, ZboundMassTol);
        syncSpeciesStorage(pair.to(),   specie, alphaSpeciesCut, ZboundMassTol);

        tmp<volScalarField> tKexp = interfacePtr->Kexp(T);
        const volScalarField& K = tKexp();

        *dmdt_[iter.key()] = K;
    }
}


template<class BasePhaseSystem>
void Foam::DiffusiveMassTransferPhaseSystem<BasePhaseSystem>::alphaTransfer
(
    SuSpTable& Su,
    SuSpTable& Sp
)
{
    return;
}


template<class BasePhaseSystem>
void Foam::DiffusiveMassTransferPhaseSystem<BasePhaseSystem>::massSpeciesTransfer
(
    const phaseModel& phase,
    volScalarField::Internal& Su,
    volScalarField::Internal& Sp,
    const word speciesName
)
{
    (void)Sp;

    auto baseSpecieName = [](const word& s) -> word
    {
        word m = IOobject::member(s);   // drops ".gas"/".liquid"

        if (m.size() && m[0] == 'X')    // tolerate old "Xxenon.gas"
        {
            m = m.substr(1);
        }

        return m;
    };

    const word target = baseSpecieName(speciesName);

    const volScalarField& alphaP = phase;
    const fvMesh& mesh = alphaP.mesh();
    const scalarField& V = mesh.V();

    const tmp<volScalarField> trhoP = phase.rho();
    const volScalarField& rhoP = trhoP();

    // Mass-transfer capacity snapshots
    // YEqns.H calls massSpeciesTransfer phase-by-phase. If the gas phase is
    // updated gas Z/Y field. These snapshots freeze receiver capacity once
    // per time index, before species source construction.
    static HashTable<scalarField> mtCapZSnapshot;
    static HashTable<scalarField> mtCapYSnapshot;
    static HashTable<label> mtCapSnapshotTime;

    static HashTable<label> capMissingWarned;

    auto snapshotKey =
    [&](const word& prefix, const word& specie, const word& phName) -> word
    {
        return word(prefix + "_" + mesh.name() + "_" + specie + "_" + phName);
    };

    auto syncAndSnapshotSpecies =
    [&]
    (
        const phaseModel& ph,
        const word& specie,
        const scalar alphaCut,
        const scalar ZboundMassTol
    ) -> void
    {
        const label timeIndex = mesh.time().timeIndex();
        const word phName = ph.name();

        const word yName0 = specie + "." + phName;
        const word yNameX = "X" + specie + "." + phName;

        const word zName0 = "Z" + specie + "." + phName;
        const word zNameX = "ZX" + specie + "." + phName;

        word yName("none");
        word zName("none");

        if (mesh.foundObject<volScalarField>(yName0))
        {
            yName = yName0;
        }
        else if (mesh.foundObject<volScalarField>(yNameX))
        {
            yName = yNameX;
        }

        if (mesh.foundObject<volScalarField>(zName0))
        {
            zName = zName0;
        }
        else if (mesh.foundObject<volScalarField>(zNameX))
        {
            zName = zNameX;
        }

        if (zName == "none" && yName == "none")
        {
            return;
        }

        const word zKey = snapshotKey("Z", specie, phName);
        const word yKey = snapshotKey("Y", specie, phName);

        const word checkKey = (zName != "none" ? zKey : yKey);

        if
        (
            mtCapSnapshotTime.found(checkKey)
         && mtCapSnapshotTime[checkKey] == timeIndex
        )
        {
            return;
        }

        const volScalarField& alphaPh = ph;

        if (zName != "none")
        {
            volScalarField& Zf =
                const_cast<volScalarField&>
                (
                    mesh.lookupObject<volScalarField>(zName)
                );

            volScalarField* YfPtr = nullptr;

            if (yName != "none")
            {
                YfPtr =
                    &const_cast<volScalarField&>
                    (
                        mesh.lookupObject<volScalarField>(yName)
                    );
            }

            const tmp<volScalarField> trhoPh = ph.rho();
            const volScalarField& rhoPh = trhoPh();

            scalarField& z = Zf.primitiveFieldRef();

            scalar dMClampLocal = 0.0;
            scalar capMassLocal = 0.0;
            scalar stockMassLocal = 0.0;

            label nHighLocal = 0;
            label nLowLocal = 0;

            forAll(z, c)
            {
                const scalar a = max(alphaPh[c], scalar(0));
                const scalar zOld = z[c];
                const scalar zNew = min(max(zOld, scalar(0)), a);

                if (zOld > a + 1e-12)
                {
                    ++nHighLocal;
                }

                if (zOld < -1e-12)
                {
                    ++nLowLocal;
                }

                if (zNew != zOld)
                {
                    dMClampLocal += rhoPh[c]*V[c]*(zNew - zOld);
                    z[c] = zNew;
                }
            }

            scalar dMClamp = dMClampLocal;
            label nHigh = nHighLocal;
            label nLow = nLowLocal;

            reduce(dMClamp, sumOp<scalar>());
            reduce(nHigh, sumOp<label>());
            reduce(nLow, sumOp<label>());

            // Redistribute any clamped mass: an absolute 1e-12 Z threshold
            // missed small undershoots and created mass at low Y
            if (mag(dMClamp) > ZboundMassTol)
            {
                if (dMClamp < 0.0)
                {
                    forAll(z, c)
                    {
                        if (alphaPh[c] > alphaCut)
                        {
                            const scalar cap = max(alphaPh[c] - z[c], scalar(0));
                            capMassLocal += rhoPh[c]*V[c]*cap;
                        }
                    }

                    scalar capMass = capMassLocal;
                    reduce(capMass, sumOp<scalar>());

                    if (capMass > ZboundMassTol)
                    {
                        const scalar frac = min((-dMClamp)/capMass, scalar(1));

                        forAll(z, c)
                        {
                            if (alphaPh[c] > alphaCut)
                            {
                                const scalar cap = max(alphaPh[c] - z[c], scalar(0));
                                const scalar dz = frac*cap;

                                z[c] += dz;
                            }
                        }
                    }
                }
                else
                {
                    forAll(z, c)
                    {
                        if (alphaPh[c] > alphaCut)
                        {
                            stockMassLocal += rhoPh[c]*V[c]*max(z[c], scalar(0));
                        }
                    }

                    scalar stockMass = stockMassLocal;
                    reduce(stockMass, sumOp<scalar>());

                    if (stockMass > ZboundMassTol)
                    {
                        const scalar frac = min(dMClamp/stockMass, scalar(1));

                        forAll(z, c)
                        {
                            if (alphaPh[c] > alphaCut)
                            {
                                const scalar dz = frac*max(z[c], scalar(0));

                                z[c] -= dz;
                            }
                        }
                    }
                }
            }

            if (YfPtr)
            {
                volScalarField& Yf = *YfPtr;
                scalarField& y = Yf.primitiveFieldRef();

                forAll(y, c)
                {
                    if (alphaPh[c] <= alphaCut)
                    {
                        y[c] = 0.0;
                    }
                    else
                    {
                        scalar yy = z[c]/max(alphaPh[c], SMALL);

                        if (yy != yy)
                        {
                            yy = 0.0;
                        }

                        y[c] = min(max(yy, scalar(0)), scalar(1));
                    }
                }

                Yf.correctBoundaryConditions();

                forAll(Zf.boundaryField(), patchi)
                {
                    fvPatchScalarField& Zp = Zf.boundaryFieldRef()[patchi];
                    const fvPatchScalarField& ap = alphaPh.boundaryField()[patchi];
                    const fvPatchScalarField& Yp = Yf.boundaryField()[patchi];

                    forAll(Zp, facei)
                    {
                        Zp[facei] = ap[facei]*Yp[facei];
                    }
                }

                Zf.correctBoundaryConditions();

                mtCapYSnapshot.set(yKey, Yf.primitiveField());
                mtCapSnapshotTime.set(yKey, timeIndex);
            }

            mtCapZSnapshot.set(zKey, Zf.primitiveField());
            mtCapSnapshotTime.set(zKey, timeIndex);
        }
        else if (yName != "none")
        {
            const volScalarField& Yf = mesh.lookupObject<volScalarField>(yName);

            mtCapYSnapshot.set(yKey, Yf.primitiveField());
            mtCapSnapshotTime.set(yKey, timeIndex);
        }
    };

    label matchedModels = 0;
    label activeCells   = 0;

    scalar mdotTargetLocal = 0.0;
    scalar mdotFromSuLocal = 0.0;

    bool dbgAny = false;

    word lastAlphaSpeciesMask("none");
    word lastAlphaControlsName("none");

    forAllConstIters(massTransferModels_, iter)
    {
        const word modelSpecie = baseSpecieName(iter()->transferSpecie());
        if (modelSpecie != target) continue;

        const phasePair& pair = iter()->pair();

        const volScalarField& alphaFrom = pair.from();
        const volScalarField& alphaTo   = pair.to();

        const bool phaseIsFrom = (phase.name() == pair.from().name());
        const bool phaseIsTo   = (phase.name() == pair.to().name());

        if (!phaseIsFrom && !phaseIsTo) continue;

        const word alphaControlsName =
            IOobject::groupName("alpha", pair.to().name());

        const dictionary& alphaControls =
            mesh.solverDict(alphaControlsName);

        const scalar alphaSpeciesCut =
            max(alphaControls.getOrDefault<scalar>("alphaSpeciesCut", 1e-12), SMALL);

        const scalar ZboundMassTol =
            alphaControls.getOrDefault<scalar>("ZboundMassTol", 1e-30);

        const scalar alphaSourceMinUser =
            alphaControls.getOrDefault<scalar>("alphaSourceMin", -1.0);

        const scalar alphaSourceMin =
            (alphaSourceMinUser > 0)
          ? max(alphaSourceMinUser, alphaSpeciesCut)
          : alphaSpeciesCut;

        const word alphaSpeciesMask =
            alphaControls.getOrDefault<word>("alphaSpeciesMask", "both");

        lastAlphaSpeciesMask = alphaSpeciesMask;
        lastAlphaControlsName = alphaControlsName;

        scalar alphaLo = alphaControls.getOrDefault<scalar>("alphaLo", 0.05);
        scalar alphaHi = alphaControls.getOrDefault<scalar>("alphaHi", 0.95);

        alphaLo = max(min(alphaLo, 0.49), 0.0);
        alphaHi = min(max(alphaHi, 0.51), 1.0);

        if (alphaHi <= alphaLo)
        {
            alphaLo = 0.05;
            alphaHi = 0.95;
        }

        const Switch dbg =
            alphaControls.getOrDefault<Switch>("debugMassSpeciesTransfer", false);

        if (dbg)
        {
            dbgAny = true;
        }

        const Switch capToCell =
            alphaControls.getOrDefault<Switch>("capToCell", true);

        const scalar YminCap =
            alphaControls.getOrDefault<scalar>("YminCap", 0.0);

        const scalar YmaxCap =
            alphaControls.getOrDefault<scalar>("YmaxCap", 1.0);

        const scalar capacitySafetyFactor =
            max
            (
                min
                (
                    alphaControls.getOrDefault<scalar>
                    (
                        "capacitySafetyFactor",
                        1.0
                    ),
                    scalar(1)
                ),
                scalar(0)
            );

        const phasePairKey k12o(pair.from().name(), pair.to().name(), true);
        const phasePairKey k12u(pair.from().name(), pair.to().name(), false);
        const phasePairKey k21o(pair.to().name(),   pair.from().name(), true);
        const phasePairKey k21u(pair.to().name(),   pair.from().name(), false);

        const volScalarField* dmdtPtr = nullptr;
        scalar sign = +1.0;

        if      (dmdt_.found(k12o)) { dmdtPtr = dmdt_[k12o]; sign = +1.0; }
        else if (dmdt_.found(k12u)) { dmdtPtr = dmdt_[k12u]; sign = +1.0; }
        else if (dmdt_.found(k21o)) { dmdtPtr = dmdt_[k21o]; sign = -1.0; }
        else if (dmdt_.found(k21u)) { dmdtPtr = dmdt_[k21u]; sign = -1.0; }
        else
        {
            continue;
        }

        const volScalarField& dmdtRaw = *dmdtPtr;

        ++matchedModels;

        // Freeze current receiver capacity before any phase species equation
        syncAndSnapshotSpecies(pair.from(), target, alphaSpeciesCut, ZboundMassTol);
        syncAndSnapshotSpecies(pair.to(),   target, alphaSpeciesCut, ZboundMassTol);

        // "both" mask: face-based conservative mapping
        if (alphaSpeciesMask == "both")
        {
            const label nCells    = mesh.nCells();
            const label nIntFaces = mesh.nInternalFaces();

            const labelList& owner = mesh.owner();
            const labelList& neigh = mesh.neighbour();

            const scalar alphaFaceSplitLocal =
                alphaControls.getOrDefault<scalar>("alphaFaceSplit", 0.5);

            const scalar interfaceAreaCut =
                alphaControls.getOrDefault<scalar>("interfaceAreaCut", 0.0);

            const Switch preferCellLocalOnMixed =
                alphaControls.getOrDefault<Switch>("preferCellLocalOnMixed", false);

            const scalar alphaBandMin =
                alphaControls.getOrDefault<scalar>("alphaBandMin", 1e-6);

            const label minMixedCells =
                alphaControls.getOrDefault<label>("minMixedCells", 50);

            const volScalarField* ifacePtr = nullptr;

            if (mesh.foundObject<volScalarField>("interfaceArea"))
            {
                ifacePtr = &mesh.lookupObject<volScalarField>("interfaceArea");
            }

            label nMixedLocal = 0;

            if (preferCellLocalOnMixed)
            {
                forAll(Su, celli)
                {
                    if (ifacePtr && (*ifacePtr)[celli] <= interfaceAreaCut) continue;

                    const scalar aF = alphaFrom[celli];
                    const scalar aT = alphaTo[celli];

                    if (aF <= alphaBandMin || aF >= 1.0 - alphaBandMin) continue;
                    if (aT <= alphaSourceMin || aT >= 1.0 - alphaBandMin) continue;

                    ++nMixedLocal;
                }
            }

            label nMixed = nMixedLocal;
            reduce(nMixed, sumOp<label>());

            if (preferCellLocalOnMixed && nMixed >= minMixedCells)
            {
                forAll(Su, celli)
                {
                    if (ifacePtr && (*ifacePtr)[celli] <= interfaceAreaCut) continue;

                    const scalar aF = alphaFrom[celli];
                    const scalar aT = alphaTo[celli];

                    if (aF <= alphaBandMin || aF >= 1.0 - alphaBandMin) continue;
                    if (aT <= alphaSourceMin || aT >= 1.0 - alphaBandMin) continue;

                    const scalar dmdtF = sign*dmdtRaw[celli];

                    if (phaseIsFrom)
                    {
                        const scalar denom = rhoP[celli]*aF;
                        if (denom <= VSMALL) continue;

                        Su[celli] -= dmdtF/denom;
                        mdotTargetLocal -= dmdtF*V[celli];
                    }
                    else
                    {
                        const scalar denom = rhoP[celli]*aT;
                        if (denom <= VSMALL) continue;

                        Su[celli] += dmdtF/denom;
                        mdotTargetLocal += dmdtF*V[celli];
                    }

                    ++activeCells;
                }

                continue;
            }

            // Branch 2: face-based mapping
            labelList ifaceFrom(nIntFaces, -1);
            labelList ifaceTo  (nIntFaces, -1);
            labelList nIfaceFrom(nCells, 0);

            label ifaceFacesLocal = 0;
            boolList hasToCell(nCells, false);

            for (label facei = 0; facei < nIntFaces; ++facei)
            {
                const label own = owner[facei];
                const label nei = neigh[facei];

                const scalar aFO = alphaFrom[own];
                const scalar aFN = alphaFrom[nei];

                const bool crossesSplit =
                    ((aFO > alphaFaceSplitLocal) != (aFN > alphaFaceSplitLocal));

                const bool brackets =
                    (min(aFO, aFN) < alphaLo && max(aFO, aFN) > alphaHi);

                if (!(crossesSplit || brackets)) continue;

                if (ifacePtr && interfaceAreaCut > SMALL)
                {
                    const scalar ia = max((*ifacePtr)[own], (*ifacePtr)[nei]);
                    if (ia <= interfaceAreaCut) continue;
                }

                bool ownFrom = false;

                if (crossesSplit)
                {
                    ownFrom = (aFO > alphaFaceSplitLocal);
                }
                else
                {
                    ownFrom = (aFO > aFN);
                }

                const label fromCell = ownFrom ? own : nei;
                const label toCell   = ownFrom ? nei : own;

                if (alphaFrom[fromCell] <= alphaSpeciesCut) continue;
                if (alphaTo[toCell]     <= alphaSourceMin) continue;

                ifaceFrom[facei] = fromCell;
                ifaceTo  [facei] = toCell;

                ++nIfaceFrom[fromCell];
                ++ifaceFacesLocal;

                hasToCell[toCell] = true;
            }

            forAll(mesh.boundary(), patchi)
            {
                if (!alphaFrom.boundaryField()[patchi].coupled()) continue;

                const fvPatch& p = mesh.boundary()[patchi];

                const fvPatchScalarField& aF_P =
                    alphaFrom.boundaryField()[patchi];

                const tmp<scalarField> taF_N = aF_P.patchNeighbourField();
                const scalarField& aF_N = taF_N();

                const fvPatchScalarField& aT_P =
                    alphaTo.boundaryField()[patchi];

                const tmp<scalarField> taT_N = aT_P.patchNeighbourField();
                const scalarField& aT_N = taT_N();

                tmp<scalarField> tia_N;
                const scalarField* ia_N_ptr = nullptr;

                if (ifacePtr)
                {
                    const fvPatchScalarField& iaP =
                        ifacePtr->boundaryField()[patchi];

                    tia_N = iaP.patchNeighbourField();
                    ia_N_ptr = &tia_N();
                }

                const labelUList& faceCells = p.faceCells();

                forAll(faceCells, fi)
                {
                    const label ownCell = faceCells[fi];

                    const scalar aFO = alphaFrom[ownCell];
                    const scalar aFN = aF_N[fi];

                    const bool crossesSplit =
                        ((aFO > alphaFaceSplitLocal) != (aFN > alphaFaceSplitLocal));

                    const bool brackets =
                        (min(aFO, aFN) < alphaLo && max(aFO, aFN) > alphaHi);

                    if (!(crossesSplit || brackets)) continue;

                    if (ifacePtr && interfaceAreaCut > SMALL)
                    {
                        const scalar iaO = (*ifacePtr)[ownCell];
                        const scalar iaN = (*ia_N_ptr)[fi];

                        if (max(iaO, iaN) <= interfaceAreaCut) continue;
                    }

                    bool ownFrom = false;

                    if (crossesSplit)
                    {
                        ownFrom = (aFO > alphaFaceSplitLocal);
                    }
                    else
                    {
                        ownFrom = (aFO > aFN);
                    }

                    if (!ownFrom) continue;

                    if (alphaFrom[ownCell] <= alphaSpeciesCut) continue;

                    if (aT_N[fi] <= alphaSourceMin) continue;

                    ++nIfaceFrom[ownCell];
                    ++ifaceFacesLocal;
                }
            }

            scalarField mdotToCell(nCells, 0.0);

            for (label facei = 0; facei < nIntFaces; ++facei)
            {
                const label fromCell = ifaceFrom[facei];
                if (fromCell < 0) continue;

                const label toCell = ifaceTo[facei];

                const label nf = nIfaceFrom[fromCell];
                if (nf <= 0) continue;

                const scalar dmdtF_from = sign*dmdtRaw[fromCell];
                const scalar mdotCell   = dmdtF_from*V[fromCell];
                const scalar mdotFace   = mdotCell/scalar(nf);

                mdotToCell[toCell] += mdotFace;
            }

            volScalarField nIfaceFromFld
            (
                IOobject
                (
                    "nIfaceFromTmp",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(dimless, Zero)
            );

            volScalarField VcellFld
            (
                IOobject
                (
                    "VcellTmp",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(dimVolume, Zero)
            );

            forAll(nIfaceFromFld, c)
            {
                nIfaceFromFld[c] = scalar(nIfaceFrom[c]);
                VcellFld[c]      = V[c];
            }

            nIfaceFromFld.correctBoundaryConditions();
            VcellFld.correctBoundaryConditions();

            volScalarField dmdtTmp
            (
                IOobject
                (
                    "dmdtTmp",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(dmdtRaw.dimensions(), Zero)
            );

            forAll(dmdtTmp, c)
            {
                dmdtTmp[c] = dmdtRaw[c];
            }

            dmdtTmp.correctBoundaryConditions();

            forAll(mesh.boundary(), patchi)
            {
                if (!alphaFrom.boundaryField()[patchi].coupled()) continue;

                const fvPatch& p = mesh.boundary()[patchi];

                const fvPatchScalarField& aF_P =
                    alphaFrom.boundaryField()[patchi];

                const tmp<scalarField> taF_N = aF_P.patchNeighbourField();
                const scalarField& aF_N = taF_N();

                const fvPatchScalarField& dmdtP =
                    dmdtTmp.boundaryField()[patchi];

                const tmp<scalarField> tdmdtN = dmdtP.patchNeighbourField();
                const scalarField& dmdtN = tdmdtN();

                const fvPatchScalarField& VP =
                    VcellFld.boundaryField()[patchi];

                const tmp<scalarField> tVN = VP.patchNeighbourField();
                const scalarField& VN = tVN();

                const fvPatchScalarField& nIFP =
                    nIfaceFromFld.boundaryField()[patchi];

                const tmp<scalarField> tnIFN = nIFP.patchNeighbourField();
                const scalarField& nIFN = tnIFN();

                tmp<scalarField> tia_N;
                const scalarField* ia_N_ptr = nullptr;

                if (ifacePtr)
                {
                    const fvPatchScalarField& iaP =
                        ifacePtr->boundaryField()[patchi];

                    tia_N = iaP.patchNeighbourField();
                    ia_N_ptr = &tia_N();
                }

                const labelUList& faceCells = p.faceCells();

                forAll(faceCells, fi)
                {
                    const label ownCell = faceCells[fi];

                    const scalar aFO = alphaFrom[ownCell];
                    const scalar aFN = aF_N[fi];

                    const bool crossesSplit =
                        ((aFO > alphaFaceSplitLocal) != (aFN > alphaFaceSplitLocal));

                    const bool brackets =
                        (min(aFO, aFN) < alphaLo && max(aFO, aFN) > alphaHi);

                    if (!(crossesSplit || brackets)) continue;

                    if (ifacePtr && interfaceAreaCut > SMALL)
                    {
                        const scalar iaO = (*ifacePtr)[ownCell];
                        const scalar iaN = (*ia_N_ptr)[fi];

                        if (max(iaO, iaN) <= interfaceAreaCut) continue;
                    }

                    bool ownFrom = false;

                    if (crossesSplit)
                    {
                        ownFrom = (aFO > alphaFaceSplitLocal);
                    }
                    else
                    {
                        ownFrom = (aFO > aFN);
                    }

                    if (ownFrom) continue;

                    if (alphaTo[ownCell] <= alphaSourceMin) continue;
                    if (aFN <= alphaSpeciesCut) continue;

                    const label nf = label(nIFN[fi]);
                    if (nf <= 0) continue;

                    const scalar dmdtF_from_nei = sign*dmdtN[fi];
                    const scalar mdotCell_nei   = dmdtF_from_nei*VN[fi];
                    const scalar mdotFace       = mdotCell_nei/scalar(nf);

                    mdotToCell[ownCell] += mdotFace;
                    hasToCell[ownCell] = true;
                }
            }

            // Capacity limiter on TO cell
            scalarField mdotToCellCapped(mdotToCell);
            scalarField scaleTo(nCells, 1.0);

            bool capAvailable = false;

            const word toYName0 = target + "." + pair.to().name();
            const word toYNameX = "X" + target + "." + pair.to().name();

            const word toZName0 = "Z" + target + "." + pair.to().name();
            const word toZNameX = "ZX" + target + "." + pair.to().name();

            const volScalarField* YtoPtr = nullptr;
            word toYNameUsed("none");

            if (mesh.foundObject<volScalarField>(toYName0))
            {
                YtoPtr = &mesh.lookupObject<volScalarField>(toYName0);
                toYNameUsed = toYName0;
            }
            else if (mesh.foundObject<volScalarField>(toYNameX))
            {
                YtoPtr = &mesh.lookupObject<volScalarField>(toYNameX);
                toYNameUsed = toYNameX;
            }

            const volScalarField* ZtoPtr = nullptr;
            word toZNameUsed("none");

            if (mesh.foundObject<volScalarField>(toZName0))
            {
                ZtoPtr = &mesh.lookupObject<volScalarField>(toZName0);
                toZNameUsed = toZName0;
            }
            else if (mesh.foundObject<volScalarField>(toZNameX))
            {
                ZtoPtr = &mesh.lookupObject<volScalarField>(toZNameX);
                toZNameUsed = toZNameX;
            }

            if (capToCell && (ZtoPtr || YtoPtr))
            {
                capAvailable = true;

                const tmp<volScalarField> trhoTo = pair.to().rho();
                const volScalarField& rhoTo = trhoTo();

                const scalar dt = max(mesh.time().deltaTValue(), SMALL);

                const scalarField* ZtoCapPtr = nullptr;
                const scalarField* YtoCapPtr = nullptr;

                word capMode("none");

                const word zSnapKey = snapshotKey("Z", target, pair.to().name());
                const word ySnapKey = snapshotKey("Y", target, pair.to().name());

                const label timeIndex = mesh.time().timeIndex();

                if
                (
                    ZtoPtr
                 && mtCapZSnapshot.found(zSnapKey)
                 && mtCapSnapshotTime.found(zSnapKey)
                 && mtCapSnapshotTime[zSnapKey] == timeIndex
                )
                {
                    ZtoCapPtr = &mtCapZSnapshot[zSnapKey];
                    capMode = "Zsnapshot";
                }
                else if
                (
                    YtoPtr
                 && mtCapYSnapshot.found(ySnapKey)
                 && mtCapSnapshotTime.found(ySnapKey)
                 && mtCapSnapshotTime[ySnapKey] == timeIndex
                )
                {
                    YtoCapPtr = &mtCapYSnapshot[ySnapKey];
                    capMode = "Ysnapshot";
                }
                else if (ZtoPtr)
                {
                    ZtoCapPtr = &(ZtoPtr->primitiveField());
                    capMode = "ZnowFallback";
                }
                else if (YtoPtr)
                {
                    YtoCapPtr = &(YtoPtr->primitiveField());
                    capMode = "YnowFallback";
                }

                label nInputLocal = 0;
                label nActiveLocal = 0;
                label nZeroedNoReceiverLocal = 0;
                label nCappedLocal = 0;

                scalar rawSignedLocal = 0.0;
                scalar cappedSignedLocal = 0.0;

                scalar rawAbsLocal = 0.0;
                scalar cappedAbsLocal = 0.0;
                scalar rejectedAbsLocal = 0.0;

                scalar minScaleLocal = GREAT;
                label worstCellLocal = -1;

                forAll(mdotToCell, c)
                {
                    const scalar mdotIn = mdotToCell[c];

                    if (mag(mdotIn) <= VSMALL)
                    {
                        scaleTo[c] = 1.0;
                        mdotToCellCapped[c] = 0.0;
                        continue;
                    }

                    ++nInputLocal;

                    rawSignedLocal += mdotIn;
                    rawAbsLocal += mag(mdotIn);

                    const scalar aTnew = max(alphaTo[c], scalar(0));

                    // The receiver must exist at the current alpha state.
                    if (aTnew <= alphaSourceMin)
                    {
                        scaleTo[c] = 0.0;
                        mdotToCellCapped[c] = 0.0;

                        ++nZeroedNoReceiverLocal;
                        ++nCappedLocal;

                        rejectedAbsLocal += mag(mdotIn);

                        minScaleLocal = 0.0;
                        worstCellLocal = c;

                        continue;
                    }

                    ++nActiveLocal;

                    const scalar aCap = aTnew;

                    scalar roomPlusZ = 0.0;
                    scalar roomMinusZ = 0.0;

                    if (ZtoCapPtr)
                    {
                        const scalar zRaw = (*ZtoCapPtr)[c];

                        const scalar zMax = aCap*YmaxCap;
                        const scalar zMin = aCap*YminCap;

                        const scalar zUsed =
                            min
                            (
                                max(zRaw, zMin),
                                zMax
                            );

                        const scalar roomPlusZRaw =
                            max(zMax - zUsed, scalar(0));

                        const scalar roomMinusZRaw =
                            max(zUsed - zMin, scalar(0));

                        roomPlusZ  = capacitySafetyFactor*roomPlusZRaw;
                        roomMinusZ = capacitySafetyFactor*roomMinusZRaw;
                    }
                    else
                    {
                        const scalar Y =
                            max
                            (
                                min((*YtoCapPtr)[c], YmaxCap),
                                YminCap
                            );

                        const scalar roomPlusZRaw =
                            aCap*max(YmaxCap - Y, scalar(0));

                        const scalar roomMinusZRaw =
                            aCap*max(Y - YminCap, scalar(0));

                        roomPlusZ  = capacitySafetyFactor*roomPlusZRaw;
                        roomMinusZ = capacitySafetyFactor*roomMinusZRaw;
                    }

                    const scalar mCapPlus = rhoTo[c]*V[c]*roomPlusZ;
                    const scalar mdotMaxPlus = mCapPlus/max(dt, SMALL);

                    const scalar mCapMinus = rhoTo[c]*V[c]*roomMinusZ;
                    const scalar mdotMinMinus = -mCapMinus/max(dt, SMALL);

                    scalar mdotC = mdotIn;

                    if (mdotIn > mdotMaxPlus)
                    {
                        mdotC = mdotMaxPlus;
                    }

                    if (mdotIn < mdotMinMinus)
                    {
                        mdotC = mdotMinMinus;
                    }

                    mdotToCellCapped[c] = mdotC;
                    scaleTo[c] = mdotC/mdotIn;

                    cappedSignedLocal += mdotC;
                    cappedAbsLocal += mag(mdotC);
                    rejectedAbsLocal += mag(mdotIn - mdotC);

                    if (mag(mdotIn - mdotC) > VSMALL)
                    {
                        ++nCappedLocal;

                        if (scaleTo[c] < minScaleLocal)
                        {
                            minScaleLocal = scaleTo[c];
                            worstCellLocal = c;
                        }
                    }
                }

                label nInput = nInputLocal;
                label nActive = nActiveLocal;
                label nZeroedNoReceiver = nZeroedNoReceiverLocal;
                label nCapped = nCappedLocal;

                scalar rawSigned = rawSignedLocal;
                scalar cappedSigned = cappedSignedLocal;

                scalar rawAbs = rawAbsLocal;
                scalar cappedAbs = cappedAbsLocal;
                scalar rejectedAbs = rejectedAbsLocal;

                scalar minScale = minScaleLocal;

                reduce(nInput, sumOp<label>());
                reduce(nActive, sumOp<label>());
                reduce(nZeroedNoReceiver, sumOp<label>());
                reduce(nCapped, sumOp<label>());

                reduce(rawSigned, sumOp<scalar>());
                reduce(cappedSigned, sumOp<scalar>());

                reduce(rawAbs, sumOp<scalar>());
                reduce(cappedAbs, sumOp<scalar>());
                reduce(rejectedAbs, sumOp<scalar>());

                reduce(minScale, minOp<scalar>());

                const bool capEvent =
                    nCapped > 0
                 || nZeroedNoReceiver > 0
                 || rejectedAbs > SMALL
                 || mag(rawAbs - cappedAbs) > SMALL;

                if (dbg && capEvent && Pstream::master())
                {
                    Info<< "MT CAP SUMMARY"
                        << " target=" << target
                        << " toPhase=" << pair.to().name()
                        << " capMode=" << capMode
                        << " Yfield=" << toYNameUsed
                        << " Zfield=" << toZNameUsed
                        << " capacitySafetyFactor=" << capacitySafetyFactor
                        << " nInput=" << nInput
                        << " nActive=" << nActive
                        << " nZeroedNoReceiver=" << nZeroedNoReceiver
                        << " nCapped=" << nCapped
                        << " rawSigned=" << rawSigned
                        << " cappedSigned=" << cappedSigned
                        << " rawAbs=" << rawAbs
                        << " cappedAbs=" << cappedAbs
                        << " rejectedAbs=" << rejectedAbs
                        << " acceptedAbsRatio=" << cappedAbs/max(rawAbs, SMALL)
                        << " rejectedAbsRatio=" << rejectedAbs/max(rawAbs, SMALL)
                        << " minScale=" << minScale
                        << nl;
                }

                if
                (
                    dbg
                 && capEvent
                 && worstCellLocal >= 0
                 && mag(minScaleLocal - minScale) <= 10*SMALL
                )
                {
                    const label c = worstCellLocal;

                    Pout<< "MT CAP WORST CELL"
                        << " target=" << target
                        << " toPhase=" << pair.to().name()
                        << " capMode=" << capMode
                        << " cell=" << c
                        << " C=" << mesh.C()[c]
                        << " V=" << V[c]
                        << " alphaTo=" << alphaTo[c]
                        << " rhoTo=" << rhoTo[c]
                        << " mdotRaw=" << mdotToCell[c]
                        << " mdotCapped=" << mdotToCellCapped[c]
                        << " scale=" << scaleTo[c];

                    if (YtoPtr)
                    {
                        Pout<< " YtoRaw=" << (*YtoPtr)[c];
                    }

                    if (ZtoPtr)
                    {
                        Pout<< " ZtoRaw=" << (*ZtoPtr)[c]
                            << " YfromZ=" << (*ZtoPtr)[c]/max(alphaTo[c], SMALL);
                    }

                    Pout<< nl;
                }
            }
            else
            {
                if (capToCell)
                {
                    const word warnKey =
                        target + "_" + pair.to().name();

                    if (!capMissingWarned.found(warnKey))
                    {
                        if (Pstream::master())
                        {
                            WarningInFunction
                                << "capToCell enabled but cannot find receiver species field "
                                << toYName0 << " or " << toYNameX
                                << ", nor conservative receiver field "
                                << toZName0 << " or " << toZNameX
                                << ". Running uncapped for this species." << nl;
                        }

                        capMissingWarned.set(warnKey, 1);
                    }
                }

                forAll(mdotToCell, c)
                {
                    scaleTo[c] = 1.0;
                    mdotToCellCapped[c] = mdotToCell[c];
                }
            }

            scalarField mdotFromAccepted(nCells, 0.0);

            volScalarField scaleToFld
            (
                IOobject
                (
                    "scaleToTmp",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(dimless, 1.0)
            );

            forAll(scaleToFld, c)
            {
                scaleToFld[c] = scaleTo[c];
            }

            scaleToFld.correctBoundaryConditions();

            for (label facei = 0; facei < nIntFaces; ++facei)
            {
                const label fromCell = ifaceFrom[facei];
                if (fromCell < 0) continue;

                const label toCell = ifaceTo[facei];

                const label nf = nIfaceFrom[fromCell];
                if (nf <= 0) continue;

                const scalar dmdtF_from = sign*dmdtRaw[fromCell];
                const scalar mdotCell   = dmdtF_from*V[fromCell];
                const scalar mdotFace   = mdotCell/scalar(nf);

                mdotFromAccepted[fromCell] += mdotFace*scaleTo[toCell];
            }

            forAll(mesh.boundary(), patchi)
            {
                if (!alphaFrom.boundaryField()[patchi].coupled()) continue;

                const fvPatch& p = mesh.boundary()[patchi];

                const fvPatchScalarField& aF_P =
                    alphaFrom.boundaryField()[patchi];

                const tmp<scalarField> taF_N = aF_P.patchNeighbourField();
                const scalarField& aF_N = taF_N();

                const fvPatchScalarField& aT_P =
                    alphaTo.boundaryField()[patchi];

                const tmp<scalarField> taT_N = aT_P.patchNeighbourField();
                const scalarField& aT_N = taT_N();

                const fvPatchScalarField& scaleP =
                    scaleToFld.boundaryField()[patchi];

                const tmp<scalarField> tScaleN = scaleP.patchNeighbourField();
                const scalarField& scaleN = tScaleN();

                tmp<scalarField> tia_N;
                const scalarField* ia_N_ptr = nullptr;

                if (ifacePtr)
                {
                    const fvPatchScalarField& iaP =
                        ifacePtr->boundaryField()[patchi];

                    tia_N = iaP.patchNeighbourField();
                    ia_N_ptr = &tia_N();
                }

                const labelUList& faceCells = p.faceCells();

                forAll(faceCells, fi)
                {
                    const label ownCell = faceCells[fi];

                    const scalar aFO = alphaFrom[ownCell];
                    const scalar aFN = aF_N[fi];

                    const bool crossesSplit =
                        ((aFO > alphaFaceSplitLocal) != (aFN > alphaFaceSplitLocal));

                    const bool brackets =
                        (min(aFO, aFN) < alphaLo && max(aFO, aFN) > alphaHi);

                    if (!(crossesSplit || brackets)) continue;

                    if (ifacePtr && interfaceAreaCut > SMALL)
                    {
                        const scalar iaO = (*ifacePtr)[ownCell];
                        const scalar iaN = (*ia_N_ptr)[fi];

                        if (max(iaO, iaN) <= interfaceAreaCut) continue;
                    }

                    bool ownFrom = false;

                    if (crossesSplit)
                    {
                        ownFrom = (aFO > alphaFaceSplitLocal);
                    }
                    else
                    {
                        ownFrom = (aFO > aFN);
                    }

                    if (!ownFrom) continue;

                    if (alphaFrom[ownCell] <= alphaSpeciesCut) continue;

                    if (aT_N[fi] <= alphaSourceMin) continue;

                    const label nf = nIfaceFrom[ownCell];
                    if (nf <= 0) continue;

                    const scalar dmdtF_from = sign*dmdtRaw[ownCell];
                    const scalar mdotCell   = dmdtF_from*V[ownCell];
                    const scalar mdotFace   = mdotCell/scalar(nf);

                    mdotFromAccepted[ownCell] += mdotFace*scaleN[fi];
                }
            }

            // donorSink drop: spread each drop's accepted donor sink over the
            // drop in proportion to its species mass (uniform fractional loss),
            // for bulk-referred closures whose depletion layer is sub-grid
            const word donorSink =
                alphaControls.getOrDefault<word>("donorSink", "interface");

            if (donorSink != "interface" && donorSink != "drop")
            {
                FatalErrorInFunction
                    << "Unsupported donorSink '" << donorSink << "' in "
                    << alphaControlsName << ". Valid choices are interface | drop"
                    << exit(FatalError);
            }

            const volScalarField* dropIdPtr =
                mesh.findObject<volScalarField>("dropId");

            const word zDonorKey = snapshotKey("Z", target, pair.from().name());

            if
            (
                phaseIsFrom
             && donorSink == "drop"
             && dropIdPtr
             && mtCapZSnapshot.found(zDonorKey)
            )
            {
                const scalarField& id = *dropIdPtr;
                const scalarField& Zd = mtCapZSnapshot[zDonorKey];

                const label nD = max(label(gMax(id)) + 1, label(0));

                scalarList Sdrop(nD, 0.0);   // accepted sink [kg/s]
                scalarList Mdrop(nD, 0.0);   // species mass [kg]

                forAll(id, c)
                {
                    const label j = label(id[c]);
                    if (j < 0) continue;

                    Sdrop[j] += mdotFromAccepted[c];

                    if (alphaFrom[c] > alphaSpeciesCut)
                    {
                        Mdrop[j] += rhoP[c]*max(Zd[c], scalar(0))*V[c];
                    }
                }

                for (label j = 0; j < nD; ++j)
                {
                    reduce(Sdrop[j], sumOp<scalar>());
                    reduce(Mdrop[j], sumOp<scalar>());
                }

                forAll(id, c)
                {
                    const label j = label(id[c]);
                    if (j < 0 || Mdrop[j] <= VSMALL) continue;

                    mdotFromAccepted[c] =
                        (alphaFrom[c] > alphaSpeciesCut)
                      ? Sdrop[j]*rhoP[c]*max(Zd[c], scalar(0))*V[c]/Mdrop[j]
                      : 0.0;
                }

                if (dbg)
                {
                    for (label j = 0; j < nD; ++j)
                    {
                        Info<< "MT DONOR DROP " << j
                            << " sink=" << Sdrop[j]
                            << " mass=" << Mdrop[j]
                            << " rate=" << Sdrop[j]/max(Mdrop[j], VSMALL)
                            << nl;
                    }
                }
            }

            const label activeBefore = activeCells;

            if (phaseIsFrom)
            {
                forAll(Su, celli)
                {
                    const scalar mdotCell = mdotFromAccepted[celli];

                    if (mag(mdotCell) <= VSMALL) continue;

                    const scalar aP = alphaFrom[celli];
                    if (aP <= alphaSpeciesCut) continue;

                    const scalar denom = rhoP[celli]*aP;
                    if (denom <= VSMALL) continue;

                    const scalar dmdtEff = mdotCell/V[celli];

                    Su[celli] -= dmdtEff/denom;
                    mdotTargetLocal -= mdotCell;

                    ++activeCells;
                }
            }
            else
            {
                forAll(Su, celli)
                {
                    const scalar mdotCell = mdotToCellCapped[celli];

                    if (mag(mdotCell) <= VSMALL) continue;

                    const scalar aP = alphaTo[celli];
                    if (aP <= alphaSourceMin) continue;

                    const scalar denom = rhoP[celli]*aP;
                    if (denom <= VSMALL) continue;

                    const scalar dmdtTo = mdotCell/V[celli];

                    Su[celli] += dmdtTo/denom;
                    mdotTargetLocal += mdotCell;

                    ++activeCells;
                }
            }

            if (dbg)
            {
                scalar sumToLocal = 0.0;
                scalar sumFromLocal = 0.0;
                scalar sumToRawLocal = 0.0;

                forAll(mdotToCell, c)
                {
                    sumToRawLocal += mdotToCell[c];
                }

                forAll(mdotToCellCapped, c)
                {
                    sumToLocal += mdotToCellCapped[c];
                }

                forAll(mdotFromAccepted, c)
                {
                    sumFromLocal += mdotFromAccepted[c];
                }

                reduce(sumToRawLocal, sumOp<scalar>());
                reduce(sumToLocal, sumOp<scalar>());
                reduce(sumFromLocal, sumOp<scalar>());

                label ifaceFacesGlobal = ifaceFacesLocal;
                reduce(ifaceFacesGlobal, sumOp<label>());

                label branchActive = activeCells - activeBefore;
                reduce(branchActive, sumOp<label>());

                const scalar mapDiff = sumToLocal - sumFromLocal;

                const scalar mapScale =
                    max
                    (
                        max(mag(sumToLocal), mag(sumFromLocal)),
                        SMALL
                    );

                const bool mapMismatch =
                    mag(mapDiff) > max(scalar(1e-30), scalar(1e-10)*mapScale);

                if (mapMismatch)
                {
                    Info<< "DEBUG both-mask: FACE-BASED mapping mismatch"
                        << " phase=" << phase.name()
                        << " species=" << speciesName
                        << " controls=" << alphaControlsName
                        << " receiver=" << pair.to().name()
                        << " donor=" << pair.from().name()
                        << " phaseIsFrom=" << phaseIsFrom
                        << " phaseIsTo=" << phaseIsTo
                        << " alphaSpeciesCut=" << alphaSpeciesCut
                        << " alphaSourceMin=" << alphaSourceMin
                        << " capToCell=" << capToCell
                        << " capAvailable=" << capAvailable
                        << " YtoField=" << toYNameUsed
                        << " ZtoField=" << toZNameUsed
                        << " ifaceFaces=" << ifaceFacesGlobal
                        << " activeCells=" << branchActive
                        << " sumToRaw=" << sumToRawLocal
                        << " sumToCapped=" << sumToLocal
                        << " sumFromAccepted=" << sumFromLocal
                        << " diff=" << mapDiff
                        << " relDiff=" << mapDiff/mapScale
                        << nl;
                }
            }

            continue;
        }

        forAll(Su, celli)
        {
            const scalar aP = alphaP[celli];

            const scalar aMin = phaseIsTo ? alphaSourceMin : alphaSpeciesCut;
            if (aP <= aMin) continue;

            if (alphaSpeciesMask == "fromOnly")
            {
                if (alphaFrom[celli] <= alphaSpeciesCut) continue;
            }

            const scalar denom = rhoP[celli]*aP;

            if (denom <= VSMALL) continue;

            const scalar dmdtF = sign*dmdtRaw[celli];

            if (phaseIsFrom)
            {
                Su[celli] -= dmdtF/denom;
                mdotTargetLocal -= dmdtF*V[celli];
            }
            else
            {
                Su[celli] += dmdtF/denom;
                mdotTargetLocal += dmdtF*V[celli];
            }

            ++activeCells;
        }
    }

    reduce(matchedModels, sumOp<label>());
    reduce(activeCells,   sumOp<label>());
    reduce(mdotTargetLocal, sumOp<scalar>());

    forAll(Su, celli)
    {
        const scalar SuCell = Su[celli];

        if (mag(SuCell) <= VSMALL) continue;

        mdotFromSuLocal += SuCell*rhoP[celli]*alphaP[celli]*V[celli];
    }

    reduce(mdotFromSuLocal, sumOp<scalar>());

    if (dbgAny && matchedModels > 0)
    {
        const scalar diff = mdotTargetLocal - mdotFromSuLocal;

        const scalar sourceScale =
            max
            (
                max(mag(mdotTargetLocal), mag(mdotFromSuLocal)),
                SMALL
            );

        const scalar rel = diff/sourceScale;

        const bool sourceMismatch =
            mag(diff) > max(scalar(1e-30), scalar(1e-10)*sourceScale);

        if (sourceMismatch)
        {
            Info<< "DEBUG massSpeciesTransfer mismatch:"
                << " phase=" << phase.name()
                << " species=" << speciesName
                << " controls=" << lastAlphaControlsName
                << " mask=" << lastAlphaSpeciesMask
                << " matchedModels=" << matchedModels
                << " activeCells=" << activeCells
                << " mdotTarget=" << mdotTargetLocal << " [kg/s]"
                << " mdotFromSu=" << mdotFromSuLocal << " [kg/s]"
                << " diff=" << diff
                << " rel=" << rel
                << nl;
        }
    }
}


template<class BasePhaseSystem>
bool Foam::DiffusiveMassTransferPhaseSystem<BasePhaseSystem>::includeVolChange()
{
    return false;
}


