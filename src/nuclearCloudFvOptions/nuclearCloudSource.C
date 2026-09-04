/*---------------------------------------------------------------------------*\
License
    This file is part of a user library, built against OpenFOAM v2606.
\*---------------------------------------------------------------------------*/

#include "nuclearCloudSource.H"
#include "fvcCurl.H"
#include "uniformDimensionedFields.H"
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
    return
        const_cast<phaseModel&>
        (
            mesh_.lookupObject<phaseSystem>("phaseProperties")
           .phases()[phaseName_]
        );
}


void Foam::fv::nuclearCloudSource::updateCarrierFields()
{
    phaseModel& p = phase();

    rhoc_ = p.rho();
    muc_ = p.thermo().mu();
    vort_ = fvc::curl(p.U());

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
                vort_,
                p.thermo().T(),
                mesh_.lookupObject<uniformDimensionedVectorField>("g")
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
        fvc::curl(phase().U())
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
