/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
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

#include "EkanayakeLiftForce.H"
#include "mathematicalConstants.H"
#include "fvcGrad.H"
#include "wallDistAddressing.H"
#include "wallPolyPatch.H"

// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::EkanayakeLiftForce<CloudType>::setWallData()
{
    const fvMesh& mesh = this->mesh();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    labelList patchIDs;
    if (patches_.empty())
    {
        DynamicList<label> ids;
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                ids.append(patchi);
            }
        }
        patchIDs.transfer(ids);
    }
    else
    {
        patchIDs = pbm.indices(patches_);
    }

    if (patchIDs.empty())
    {
        FatalErrorInFunction
            << "No wall patches for " << this->typeName
            << exit(FatalError);
    }

    // Nearest-face addressing over the selected patches (computed once,
    // the mesh is static)
    const auto& wDist = wallDistAddressing::New
    (
        "EkanayakeLiftForce:walls",
        mesh,
        patchIDs
    );

    volVectorField n
    (
        IOobject
        (
            "EkanayakeLiftForce:nWall",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedVector(dimless, Zero)
    );

    volVectorField Cf
    (
        IOobject
        (
            "EkanayakeLiftForce:CfWall",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedVector(dimLength, Zero)
    );

    for (const label patchi : wDist.patchIDs())
    {
        n.boundaryFieldRef()[patchi] == -mesh.boundary()[patchi].nf();
        Cf.boundaryFieldRef()[patchi] == mesh.boundary()[patchi].Cf();
    }

    wDist.map(n, mapDistribute::transform());
    wDist.map(Cf, mapDistribute::transformPosition());

    nWall_ = n.primitiveField();
    CfWall_ = Cf.primitiveField();

    Info<< "    " << this->typeName << ": walls "
        << UIndirectList<word>(pbm.names(), patchIDs)
        << ", dMin " << dMin_
        << ", hinderedMigration " << hinderedMigration_
        << ", ReMax " << ReMax_
        << ", slip " << (terminalSlip_ ? "terminal" : "relative") << endl;
}


template<class CloudType>
Foam::vector Foam::EkanayakeLiftForce<CloudType>::terminalSlip
(
    const scalar d,
    const scalar rhop,
    const scalar rhoc,
    const scalar muc
) const
{
    // Settling velocity relative to the fluid: buoyancy-corrected gravity
    // balanced by the Schiller-Naumann drag of sphereDrag,
    // u_t = u_Stokes/(1 + 0.15 Re^0.687), Re = rho_c u_t d/mu_c
    const vector& g = this->owner().g().value();
    const scalar magg = mag(g);

    if (magg < VSMALL)
    {
        return Zero;
    }

    const scalar uStokes = (rhop - rhoc)*magg*sqr(d)/(18.0*muc);

    scalar ut = uStokes;
    for (label i = 0; i < 50; ++i)
    {
        const scalar utNew =
            uStokes/(1.0 + 0.15*pow(rhoc*mag(ut)*d/muc, 0.687));

        if (mag(utNew - ut) <= 1e-10*mag(ut))
        {
            ut = utNew;
            break;
        }
        ut = utNew;
    }

    return ut*g/magg;
}


template<class CloudType>
Foam::scalar Foam::EkanayakeLiftForce<CloudType>::J(const scalar epsilon)
{
    // Legendre & Magnaudet (1998); J -> 2.255 (Saffman) for epsilon -> inf
    return 2.255*pow(1.0 + 0.2/sqr(max(epsilon, 1e-6)), -1.5);
}


template<class CloudType>
Foam::scalar Foam::EkanayakeLiftForce<CloudType>::CW0(const scalar Lu)
{
    // Vasseur & Cox (1977) fitted by Shi & Rzehak (2020, eq 20)
    return
        Lu < 10.0
      ? 9.0/8.0/(1.0 + 0.13*Lu*(Lu + 0.53))
      : 8.94*pow(Lu, -2.09);
}


template<class CloudType>
Foam::scalar Foam::EkanayakeLiftForce<CloudType>::CW
(
    const scalar LR,
    const scalar Rep
)
{
    // Zeng et al. (2009) with eq 20 as the Re_p << 1 base
    // (Shi & Rzehak 2020, eq 22)
    const scalar A = 1.0 + 0.329*Rep + 0.00485*sqr(Rep);
    const scalar g = -0.9*tanh(0.022*Rep);
    const scalar delta = 0.5*(LR - 1.0);
    const scalar deltac = 3.0*exp(-0.17*pow(Rep, 0.7));
    const scalar C101 = 0.313 + 0.812*exp(-0.125*pow(Rep, 0.77));

    const scalar c = CW0(0.5*LR*Rep)*pow(A, pow(0.5*LR, g));
    const scalar c101 = CW0(0.5*1.01*Rep)*pow(A, pow(0.5*1.01, g));

    return c + (C101 - c101)*exp(-11.0*pow(max(delta, 0.0)/deltac, 1.2));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::EkanayakeLiftForce<CloudType>::EkanayakeLiftForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& forceType
)
:
    LiftForce<CloudType>(owner, mesh, dict, forceType),
    dMin_(this->coeffs().template getOrDefault<scalar>("dMin", 0)),
    hinderedMigration_
    (
        this->coeffs().template getOrDefault<bool>("hinderedMigration", true)
    ),
    ReMax_(this->coeffs().template getOrDefault<scalar>("ReMax", 0.1)),
    patches_
    (
        this->coeffs().template getOrDefault<wordRes>("patches", wordRes())
    ),
    terminalSlip_(false),
    nWall_(),
    CfWall_(),
    gradUcInterpPtr_(nullptr),
    warned_(false)
{
    const word slip
    (
        this->coeffs().template getOrDefault<word>("slip", "relative")
    );

    if (slip == "terminal")
    {
        terminalSlip_ = true;
    }
    else if (slip != "relative")
    {
        FatalIOErrorInFunction(this->coeffs())
            << "slip " << slip << " unknown: relative or terminal"
            << exit(FatalIOError);
    }

    setWallData();
}


template<class CloudType>
Foam::EkanayakeLiftForce<CloudType>::EkanayakeLiftForce
(
    const EkanayakeLiftForce& lf
)
:
    LiftForce<CloudType>(lf),
    dMin_(lf.dMin_),
    hinderedMigration_(lf.hinderedMigration_),
    ReMax_(lf.ReMax_),
    patches_(lf.patches_),
    terminalSlip_(lf.terminalSlip_),
    nWall_(lf.nWall_),
    CfWall_(lf.CfWall_),
    gradUcInterpPtr_(nullptr),
    warned_(false)
{}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

template<class CloudType>
void Foam::EkanayakeLiftForce<CloudType>::cacheFields(const bool store)
{
    LiftForce<CloudType>::cacheFields(store);

    static word resultName("gradUcDt");

    volTensorField* resultPtr =
        this->mesh().template getObjectPtr<volTensorField>(resultName);

    if (store)
    {
        if (!resultPtr)
        {
            const volVectorField& Uc = this->mesh().template
                lookupObject<volVectorField>(this->UName_);

            resultPtr = new volTensorField(resultName, fvc::grad(Uc));

            resultPtr->store();
        }

        gradUcInterpPtr_.reset
        (
            interpolation<tensor>::New
            (
                this->owner().solution().interpolationSchemes(),
                *resultPtr
            ).ptr()
        );
    }
    else
    {
        gradUcInterpPtr_.clear();

        if (resultPtr)
        {
            resultPtr->checkOut();
        }
    }
}


template<class CloudType>
Foam::forceSuSp Foam::EkanayakeLiftForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    using constant::mathematical::pi;

    forceSuSp value(Zero);

    const scalar d = p.d();

    if (d < dMin_)
    {
        return value;
    }

    const scalar rhoc = td.rhoc();
    const scalar nuc = muc/rhoc;
    const scalar a = 0.5*d;

    // Unbounded shear lift (McLaughlin 1991 with Legendre & Magnaudet's J),
    // vector form of LiftForce: F = (m/rho_p) rho_c Cl (U_c - U_p) x curl U_c
    const vector curlUc =
        this->curlUcInterp().interpolate(p.coordinates(), p.currentTetIndices());

    const vector Us =
        terminalSlip_
      ? terminalSlip(d, p.rho(), rhoc, muc)
      : vector(p.U() - td.Uc());
    const scalar Rew = mag(curlUc)*sqr(d)/nuc;
    const scalar epsUb = sqrt(mag(curlUc)*nuc)/(mag(Us) + ROOTVSMALL);
    const scalar Cl = 27.0*J(epsUb)/(2.0*sqr(pi)*sqrt(Rew + ROOTVSMALL));

    const vector Fub = mass/p.rho()*rhoc*Cl*(-Us ^ curlUc);

    // Nearest wall
    const label celli = p.cell();
    const vector& n = nWall_[celli];

    if (magSqr(n) < 0.25)
    {
        // No wall reached from this cell: unbounded lift
        value.Su() = Fub;
        return value;
    }

    const scalar l = (p.position() - CfWall_[celli]) & n;
    const scalar lStar = max(l/a, 1.2);
    const scalar s = 1.0/lStar;

    // Slip along the wall and the wall shear along it
    const tensor gradUc =
        gradUcInterpPtr_().interpolate
        (
            p.coordinates(),
            p.currentTetIndices()
        );

    const vector shearVec = n & gradUc;     // dU/dn
    const vector UsPar = Us - (Us & n)*n;
    const scalar magUsPar = mag(UsPar);

    vector t(Zero);
    if (magUsPar > ROOTVSMALL)
    {
        t = UsPar/magUsPar;
    }
    else
    {
        // No slip along the wall: shear direction (C_L,1 needs |gamma| only)
        const vector shearPar = shearVec - (shearVec & n)*n;
        t = shearPar/(mag(shearPar) + ROOTVSMALL);
    }

    const scalar gamma = shearVec & t;      // dU_t/dn
    const scalar Reg = mag(gamma)*sqr(a)/nuc;
    const scalar Res = magUsPar*a/nuc;

    const scalar x = lStar*sqrt(Reg);       // l/L_G
    const scalar y = lStar*Res;             // l/L_S
    const scalar eps = sqrt(Reg)/(Res + ROOTVSMALL);

    const bool finiteRe = (Reg > ReMax_ || Res > ReMax_);

    // f: factor on the unbounded lift; C*: wall-normal coefficients
    // (debug output only for the branch not taken)
    scalar f = 0, Fn = 0;
    scalar CL1 = 0, CL2in = 0, CL3 = 0, CWL = 0, fL = 0;

    if (finiteRe)
    {
        // Takemura et al. (2009) decomposition with the correlations of
        // Shi & Rzehak (2020, eq 24, Table 3); Re_p = |u_s| d/nu
        const scalar Rep = 2.0*Res;

        if (!warned_ && Rep > 4.0)
        {
            WarningInFunction
                << "Re_p " << Rep << " above 4 (d " << d
                << "): outside the range of the finite-Re lift of "
                << this->typeName << " (warned once)" << endl;
            warned_ = true;
        }

        // Wall modification of the shear lift (Shi & Rzehak eq 29)
        const scalar magJ = mag(J(eps));
        f =
            x < 15.0
          ? 1.0
          - exp(-11.0/96.0*sqr(pi)*70.0/(70.0 + pow(x, 1.378))*x/magJ)
          : 1.0 - 1.8778*pow(x, -5.0/3.0)/magJ;

        // Wall lift (Zeng et al. 2009) and its reduction by the shear
        // (Takemura et al. 2009, Shi & Rzehak eq 30)
        CWL = CW(lStar, Rep);
        fL = exp(-0.22*pow(eps, 0.8)*pow(x, 2.5));

        Fn = CWL*fL*pi/8.0*sqr(d)*rhoc*sqr(magUsPar);
    }
    else
    {
        // Wall reduction of the shear lift (Takemura et al. 2009)
        f = 1.0 - exp(-11.0/96.0*sqr(pi)*x/J(eps));

        // Zero-slip lift, freely rotating (Ekanayake et al. 2020)
        const scalar f1 = 0.925*exp(-0.35*Reg) - 0.0135*exp(-7000.0*Reg);
        const scalar f2 = 1.0 + sqrt(Reg);
        CL1 =
            f1*1.982*exp(-0.115*sqr(x) - 0.2771*x)
          + f2*(0.89934*s - 1.9610*sqr(s) + 1.0216*pow3(s));

        // Inner slip-shear lift, freely rotating (Cherukat & McLaughlin 1995)
        CL2in = -2.6729 - 0.8373*s + 0.4683*sqr(s);

        // Slip (wall) lift, freely rotating (Ekanayake et al. 2021, eq 4.4b)
        CL3 =
            18.0*pi/(32.0 + 2.0*y + 3.8*sqr(y) + 0.049*pow3(y))
          + 0.4757*s - 1.268*sqr(s) + 0.683*pow3(s);

        Fn =
            sqr(muc)/rhoc
           *(
                CL1*sqr(Reg)
              + sign(gamma)*CL2in*Reg*Res
              + CL3*sqr(Res)
            );
    }

    vector F = f*Fub + Fn*n;

    if (hinderedMigration_)
    {
        // Wall-normal drag hindrance: 6 pi/C_D,perp (Brenner)
        const scalar hPerp =
            max
            (
                1.0 - 9.0/8.0*s + 0.5*pow3(s)
              - 135.0/256.0*pow4(s) - 0.125*pow5(s),
                0.0
            );

        F += (hPerp - 1.0)*(F & n)*n;
    }

    if (debug)
    {
        Pout<< "EkanayakeLift: d " << d << " x " << p.position()
            << " n " << n << " l " << l << " Us " << Us
            << " gamma " << gamma << " Reg " << Reg << " Res " << Res
            << " f " << f << " CL1 " << CL1 << " CL2in " << CL2in
            << " CL3 " << CL3 << " CW " << CWL << " fL " << fL
            << " Fn " << Fn << " Fub " << Fub
            << " F " << F << endl;
    }

    value.Su() = F;

    return value;
}


// ************************************************************************* //
