/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include "ForceMagnitudes.H"
#include "fvcCurl.H"
#include "fvcDdt.H"
#include "fvcGrad.H"
#include "interpolation.H"
#include "mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ForceMagnitudes<CloudType>::ForceMagnitudes
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    UName_(dict.getOrDefault<word>("U", "U")),
    rhoc_(dict.getScalar("rhoc")),
    nuc_(dict.getScalar("nuc")),
    Cvm_(dict.getOrDefault<scalar>("Cvm", 0.5))
{}


template<class CloudType>
Foam::ForceMagnitudes<CloudType>::ForceMagnitudes
(
    const ForceMagnitudes<CloudType>& fm
)
:
    CloudFunctionObject<CloudType>(fm),
    UName_(fm.UName_),
    rhoc_(fm.rhoc_),
    nuc_(fm.nuc_),
    Cvm_(fm.Cvm_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ForceMagnitudes<CloudType>::postEvolve
(
    const typename parcelType::trackingData& td
)
{
    auto& c = this->owner();

    const bool haveParticles = c.size();

    if (!returnReduceOr(haveParticles))
    {
        return;
    }

    const fvMesh& mesh = c.mesh();
    const volVectorField& Uc = mesh.template lookupObject<volVectorField>(UName_);

    // Same auxiliary fields as LiftForce/PressureGradientForce, computed
    // locally here so this function does not depend on those forces being
    // active in "particleForces". Named to match the LiftForce/
    // PressureGradientForce field names ("curlUcDt"/"DUcDt") so they pick
    // up the same interpolationSchemes entries as those force classes.
    const volVectorField curlUc("curlUcDt", fvc::curl(Uc));
    const volVectorField DUcDt("DUcDt", fvc::ddt(Uc) + (Uc & fvc::grad(Uc)));

    autoPtr<interpolation<vector>> UcInterp
    (
        interpolation<vector>::New(c.solution().interpolationSchemes(), Uc)
    );
    autoPtr<interpolation<vector>> curlUcInterp
    (
        interpolation<vector>::New
        (
            c.solution().interpolationSchemes(),
            curlUc
        )
    );
    autoPtr<interpolation<vector>> DUcDtInterp
    (
        interpolation<vector>::New
        (
            c.solution().interpolationSchemes(),
            DUcDt
        )
    );

    const scalar muc = rhoc_*nuc_;
    const vector g = c.g().value();

    auto initField = [&](const word& name) -> IOField<scalar>&
    {
        auto* ptr = c.template getObjectPtr<IOField<scalar>>(name);

        if (!ptr)
        {
            ptr = new IOField<scalar>
            (
                IOobject
                (
                    name,
                    c.time().timeName(),
                    c,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::REGISTER
                )
            );

            ptr->store();
        }

        ptr->resize(c.size());

        return *ptr;
    };

    auto& FDrag = initField("FDrag");
    auto& FGrav = initField("FGrav");
    auto& FLift = initField("FLift");
    auto& FPg   = initField("FPg");
    auto& FVm   = initField("FVm");

    label parceli = 0;

    for (const parcelType& p : c)
    {
        const vector Uc_p =
            UcInterp().interpolate(p.coordinates(), p.currentTetIndices());
        const vector curlUc_p =
            curlUcInterp().interpolate(p.coordinates(), p.currentTetIndices());
        const vector DUcDt_p =
            DUcDtInterp().interpolate(p.coordinates(), p.currentTetIndices());

        const vector Urel = Uc_p - p.U();
        const scalar d = p.d();
        const scalar rhop = p.rho();
        const scalar mass = p.mass();

        const scalar Re = mag(Urel)*d/nuc_;

        // Drag - Foam::SphereDragForce (AOB:Eq. 34/35)
        {
            const scalar CdRe =
                (Re > 1000)
              ? 0.424*Re
              : 24.0*(1.0 + (1.0/6.0)*pow(Re, 2.0/3.0));

            const scalar Sp = mass*0.75*muc*CdRe/(rhop*sqr(d));

            FDrag[parceli] = mag(Sp*Urel);
        }

        // Gravity + buoyancy - Foam::GravityForce
        FGrav[parceli] = mag(mass*g*(1.0 - rhoc_/rhop));

        // Saffman-Mei lift - Foam::SaffmanMeiLiftForce / Foam::LiftForce
        {
            const scalar Rew =
                rhoc_*mag(curlUc_p)*sqr(d)/(muc + ROOTVSMALL);
            const scalar beta = 0.5*(Rew/(Re + ROOTVSMALL));
            const scalar alpha = 0.3314*sqrt(beta);
            const scalar f = (1.0 - alpha)*exp(-0.1*Re) + alpha;

            const scalar Cld =
                (Re < 40) ? 6.46*f : 6.46*0.0524*sqrt(beta*Re);

            const scalar Cl =
                3.0/(mathematical::twoPi*sqrt(Rew + ROOTVSMALL))*Cld;

            FLift[parceli] = mag(mass/rhop*rhoc_*Cl*(Urel ^ curlUc_p));
        }

        // Pressure gradient - Foam::PressureGradientForce
        // Virtual (added) mass - Foam::VirtualMassForce (= Cvm * pg force)
        {
            const vector FpgVec = mass*rhoc_/rhop*DUcDt_p;

            FPg[parceli] = mag(FpgVec);
            FVm[parceli] = mag(Cvm_*FpgVec);
        }

        ++parceli;
    }

    if (c.time().writeTime())
    {
        FDrag.write(haveParticles);
        FGrav.write(haveParticles);
        FLift.write(haveParticles);
        FPg.write(haveParticles);
        FVm.write(haveParticles);
    }
}


// ************************************************************************* //
