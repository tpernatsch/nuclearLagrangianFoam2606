/*---------------------------------------------------------------------------*\
License
    This file is part of a user library, built against OpenFOAM v2606.
\*---------------------------------------------------------------------------*/

#include "nuclearCloudSource.H"
#include "fvcCurl.H"
#include "gravityMeshObject.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(nuclearCloudSource, 0);

    addToRunTimeSelectionTable
    (
        option,
        nuclearCloudSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

Foam::phaseModel& Foam::fv::nuclearCloudSource::phase() const
{
    // Avoids linking reactingMultiphaseSystem for phaseSystem's RTTI -
    // look it up as the always-linked IOdictionary base and cast down
    const IOdictionary& phaseProperties =
        mesh_.lookupObject<IOdictionary>("phaseProperties");

    const phaseSystem& fluid =
        static_cast<const phaseSystem&>(phaseProperties);

    return const_cast<phaseModel&>(fluid.phases()[phaseName_]);
}


void Foam::fv::nuclearCloudSource::updateCarrierFields()
{
    phaseModel& p = phase();

    rhoc_ = p.rho();
    muc_ = p.thermo().mu();

    if (!cloud_)
    {
        cloud_.reset
        (
            new basicNuclearCloud
            (
                cloudName_,
                rhoc_,
                p.URef(),
                muc_,
                p.thermo().T(),
                meshObjects::gravity::New(mesh_.time())
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::fv::nuclearCloudSource::nuclearCloudSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    phaseName_(coeffs_.lookup("phase")),
    cloudName_(coeffs_.getOrDefault<word>("cloudName", name_ + "Cloud")),
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
        phase().rho()
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
        phase().thermo().mu()
    ),
    cloud_(nullptr)
{
    fieldNames_.setSize(1, IOobject::groupName("U", phaseName_));
    resetApplied();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::nuclearCloudSource::correct(volVectorField& U)
{
    // Skip the very first call: some carrier fields (e.g. alphaPhi-derived
    // quantities) are not meaningfully evaluable until after the first
    // phase/pressure solve has run once, matching the same guard
    // flotationFoam's fv::flotation::correct() uses.
    if (mesh_.time().timeIndex() <= 1)
    {
        return;
    }

    updateCarrierFields();

    cloud_->storeGlobalPositions();
    cloud_->evolve();
}


bool Foam::fv::nuclearCloudSource::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        return true;
    }

    return false;
}


// ************************************************************************* //
