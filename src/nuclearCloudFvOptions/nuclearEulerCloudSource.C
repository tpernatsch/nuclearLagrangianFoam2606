/*---------------------------------------------------------------------------*\
License
    This file is part of a user library, built against OpenFOAM v2606.
\*---------------------------------------------------------------------------*/

#include "nuclearEulerCloudSource.H"
#include "phaseCompressibleTurbulenceModel.H"
#include "fvcCurl.H"
#include "gravityMeshObject.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(nuclearEulerCloudSource, 0);

    addToRunTimeSelectionTable
    (
        option,
        nuclearEulerCloudSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

Foam::phaseModel& Foam::fv::nuclearEulerCloudSource::phaseByName
(
    const word& phaseName
) const
{
    // Avoids linking reactingMultiphaseSystem for phaseSystem's RTTI -
    // look it up as the always-linked IOdictionary base and cast down
    const IOdictionary& phaseProperties =
        mesh_.lookupObject<IOdictionary>("phaseProperties");

    const phaseSystem& fluid =
        static_cast<const phaseSystem&>(phaseProperties);

    return const_cast<phaseModel&>(fluid.phases()[phaseName]);
}


Foam::phaseModel& Foam::fv::nuclearEulerCloudSource::phase() const
{
    return phaseByName(phaseName_);
}


Foam::phaseModel& Foam::fv::nuclearEulerCloudSource::phase2() const
{
    return phaseByName(phase2Name_);
}


const Foam::phaseCompressibleTurbulenceModel&
Foam::fv::nuclearEulerCloudSource::turbulenceByName
(
    const word& phaseName
) const
{
    return mesh_.lookupObject<phaseCompressibleTurbulenceModel>
    (
        IOobject::groupName(turbulenceModel::propertiesName, phaseName)
    );
}


bool Foam::fv::nuclearEulerCloudSource::isRASByName
(
    const word& phaseName
) const
{
    return
        turbulenceByName(phaseName).get<word>("simulationType")
     == "RAS";
}


void Foam::fv::nuclearEulerCloudSource::updateCarrierFields()
{
    phaseModel& p = phase();
    phaseModel& p2 = phase2();

    // phaseModel derives from volScalarField (it *is* the alpha field);
    // phaseModel::alpha() is a same-named but unrelated thermo interface
    // method (thermal diffusivity, kg/m/s), not the volume fraction.
    alpha1c_ = p;
    rhoc_ = p.rho();
    muc_ = p.thermo().mu();
    k1c_ = turbulenceByName(phaseName_).k();
    epsilon1c_ = turbulenceByName(phaseName_).epsilon();

    alpha2c_ = p2;
    rho2c_ = p2.rho();
    mu2c_ = p2.thermo().mu();
    k2c_ = turbulenceByName(phase2Name_).k();
    epsilon2c_ = turbulenceByName(phase2Name_).epsilon();
    d2c_ = p2.d();

    if (!cloud_)
    {
        cloud_.reset
        (
            new basicNuclearEulerCloud
            (
                cloudName_,
                rhoc_,
                p.URef(),
                muc_,
                p.thermo().T(),
                meshObjects::gravity::New(mesh_.time()),
                alpha1c_,
                k1c_,
                epsilon1c_,
                isRASByName(phaseName_),
                eulerCarrierPhaseFields
                (
                    alpha2c_,
                    p2.URef(),
                    rho2c_,
                    mu2c_,
                    p2.thermo().T(),
                    k2c_,
                    epsilon2c_,
                    d2c_,
                    isRASByName(phase2Name_)
                )
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::fv::nuclearEulerCloudSource::nuclearEulerCloudSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    phaseName_(coeffs_.lookup("phase")),
    phase2Name_(coeffs_.lookup("phase2")),
    cloudName_(coeffs_.getOrDefault<word>("cloudName", name_ + "Cloud")),
    alpha1c_
    (
        IOobject
        (
            name_ + ":alpha1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phaseByName(phaseName_)
    ),
    rhoc_
    (
        IOobject
        (
            name_ + ":rho",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phaseByName(phaseName_).rho()
    ),
    muc_
    (
        IOobject
        (
            name_ + ":mu",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phaseByName(phaseName_).thermo().mu()
    ),
    k1c_
    (
        IOobject
        (
            name_ + ":k1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        turbulenceByName(phaseName_).k()
    ),
    epsilon1c_
    (
        IOobject
        (
            name_ + ":epsilon1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        turbulenceByName(phaseName_).epsilon()
    ),
    alpha2c_
    (
        IOobject
        (
            name_ + ":alpha2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phaseByName(phase2Name_)
    ),
    rho2c_
    (
        IOobject
        (
            name_ + ":rho2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phaseByName(phase2Name_).rho()
    ),
    mu2c_
    (
        IOobject
        (
            name_ + ":mu2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phaseByName(phase2Name_).thermo().mu()
    ),
    k2c_
    (
        IOobject
        (
            name_ + ":k2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        turbulenceByName(phase2Name_).k()
    ),
    epsilon2c_
    (
        IOobject
        (
            name_ + ":epsilon2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        turbulenceByName(phase2Name_).epsilon()
    ),
    d2c_
    (
        IOobject
        (
            name_ + ":d2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phaseByName(phase2Name_).d()
    ),
    cloud_(nullptr)
{
    fieldNames_.setSize(1, IOobject::groupName("U", phaseName_));
    resetApplied();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::nuclearEulerCloudSource::correct(volVectorField& U)
{
    // Skip the very first call, carrier fields not yet meaningful
    if (mesh_.time().timeIndex() <= 1)
    {
        return;
    }

    updateCarrierFields();

    cloud_->storeGlobalPositions();
    cloud_->evolve();
}


bool Foam::fv::nuclearEulerCloudSource::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        return true;
    }

    return false;
}


// ************************************************************************* //
