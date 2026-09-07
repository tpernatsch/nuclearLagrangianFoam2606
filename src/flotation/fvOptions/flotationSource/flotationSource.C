/*---------------------------------------------------------------------------*\
License
    This file is part of a user library, built against OpenFOAM v2606.
\*---------------------------------------------------------------------------*/

#include "flotationSource.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(flotationSource, 0);

    addToRunTimeSelectionTable
    (
        option,
        flotationSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::flotationSource::flotationSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh)
{
    const phaseSystem& fluid =
        mesh.lookupObject<phaseSystem>("phaseProperties");

    if (coeffs_.found("dict"))
    {
        system_.reset
        (
            new flotationSystem(fluid, mesh, word(coeffs_.lookup("dict")))
        );
    }
    else
    {
        system_.reset(new flotationSystem(fluid, mesh));
    }

    // Piggyback on the first interface phase's U correct() call as a
    // once-per-outer-iteration clock tick (mirrors fv::flotation::correct()
    // in the original Foundation fvModel, which is called unconditionally
    // once per iteration - ESI's fv::option has no such field-agnostic
    // hook, only per-field correct(field), so we tie to one specific,
    // definitely-moving field instead).
    fieldNames_.setSize
    (
        1,
        IOobject::groupName("U", system_->interface().phase1().name())
    );
    resetApplied();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::flotationSource::correct(volVectorField& U)
{
    if (mesh_.time().timeIndex() > 1)
    {
        system_->solve();
    }
}


bool Foam::fv::flotationSource::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        return true;
    }

    return false;
}


// ************************************************************************* //
