/*---------------------------------------------------------------------------*\
License
    This file is part of a user library, built against OpenFOAM v2606.
\*---------------------------------------------------------------------------*/

#include "nuclearCloudSourceSingle.H"
#include "fvcCurl.H"
#include "gravityMeshObject.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(nuclearCloudSourceSingle, 0);

    addToRunTimeSelectionTable
    (
        option,
        nuclearCloudSourceSingle,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::fv::nuclearCloudSourceSingle::nuclearCloudSourceSingle
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    cloudName_(coeffs_.getOrDefault<word>("cloudName", name_ + "Cloud")),
    rho0_("rho", dimDensity, coeffs_.get<scalar>("rho")),
    mu0_
    (
        "mu",
        dimDynamicViscosity,
        coeffs_.get<scalar>("rho")*coeffs_.get<scalar>("nu")
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
        mesh_,
        rho0_
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
        mesh_,
        mu0_
    ),
    cloud_(nullptr),
    lastTimeIndex_(-1)
{
    fieldNames_.setSize(1, "U");
    resetApplied();

    // Build the cloud now, at the start time: it then reads the parcels and
    // its saved state (uniform/lagrangian/<cloud>, e.g. what the injection
    // models have already released) of a restart. Built later, in correct(),
    // the time has already advanced: the cloud would read nothing, start
    // empty, and the injection would release everything again. U and T exist
    // here: the solvers create the fvOptions after their fields.
    cloud_.reset
    (
        new basicNuclearCloud
        (
            cloudName_,
            rhoc_,
            mesh_.lookupObject<volVectorField>("U"),
            muc_,
            mesh_.lookupObject<volScalarField>("T"),
            meshObjects::gravity::New(mesh_.time())
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::nuclearCloudSourceSingle::correct(volVectorField& U)
{
    // Evolve once per time step, on the first call (after the momentum
    // predictor); later calls from the pressure correctors are ignored.
    // No "settle first" skip of the first step: the cloud exists from the
    // start (and a restart must not lose a step of particle motion).
    if (mesh_.time().timeIndex() == lastTimeIndex_)
    {
        return;
    }
    lastTimeIndex_ = mesh_.time().timeIndex();

    cloud_->storeGlobalPositions();
    cloud_->evolve();
}


bool Foam::fv::nuclearCloudSourceSingle::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        return true;
    }

    return false;
}


// ************************************************************************* //
