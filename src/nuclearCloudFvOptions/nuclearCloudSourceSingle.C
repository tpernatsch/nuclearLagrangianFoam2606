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

void Foam::fv::nuclearCloudSourceSingle::correct(volVectorField& U)
{
    // Same "settle first" guard as nuclearCloudSource/kinematicCloudSource
    if (mesh_.time().timeIndex() <= 1)
    {
        return;
    }

    vort_ = fvc::curl(U);

    if (!cloud_)
    {
        cloud_.reset
        (
            new basicNuclearCloud
            (
                cloudName_,
                rhoc_,
                U,
                muc_,
                vort_,
                mesh_.lookupObject<volScalarField>("T"),
                meshObjects::gravity::New(mesh_.time())
            )
        );
    }

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
