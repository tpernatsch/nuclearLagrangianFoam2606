/*---------------------------------------------------------------------------*\
License
    This file is part of a user library, built against OpenFOAM v2606.
\*---------------------------------------------------------------------------*/

#include "kinematicCloudSource.H"
#include "fvcCurl.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(kinematicCloudSource, 0);

    addToRunTimeSelectionTable
    (
        option,
        kinematicCloudSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::fv::kinematicCloudSource::kinematicCloudSource
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
    g_
    (
        "g",
        dimAcceleration,
        coeffs_.getOrDefault<vector>("g", Zero)
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
    vort_
    (
        IOobject
        (
            name_ + ":vort",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::curl(mesh_.lookupObject<volVectorField>("U"))
    ),
    cloud_(nullptr)
{
    fieldNames_.setSize(1, "U");
    resetApplied();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::kinematicCloudSource::correct(volVectorField& U)
{
    // Skip the very first call, matching fv::flotation::correct() and
    // nuclearCloudSource::correct() - keeps the same "settle first" guard
    // consistent across all three, though here it mainly just avoids
    // injecting on an as-yet-unsolved flow field.
    if (mesh_.time().timeIndex() <= 1)
    {
        return;
    }

    vort_ = fvc::curl(U);

    if (!cloud_)
    {
        cloud_.reset
        (
            new basicKinematicCloud
            (
                cloudName_,
                rhoc_,
                U,
                muc_,
                vort_,
                g_
            )
        );
    }

    cloud_->storeGlobalPositions();
    cloud_->evolve();
}


bool Foam::fv::kinematicCloudSource::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        return true;
    }

    return false;
}


// ************************************************************************* //
