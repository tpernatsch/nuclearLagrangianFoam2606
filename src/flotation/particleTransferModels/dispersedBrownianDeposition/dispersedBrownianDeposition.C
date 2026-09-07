#include "dispersedBrownianDeposition.H"
#include "particleModel.H"
#include "flotationSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace particleTransferModels
{
    defineTypeNameAndDebug(dispersedBrownianDeposition, 0);
    addToRunTimeSelectionTable
    (
        particleTransferModel,
        dispersedBrownianDeposition,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particleTransferModels::dispersedBrownianDeposition::
dispersedBrownianDeposition
(
    const particleModel& model,
    const dictionary& dict
)
:
    dispersedParticleTransferModel(model, dict),
    Sh_
    (
        SherwoodNumber::New
        (
            model.system().mesh(),
            word(dict.lookup("Sh"))
        ).ptr()
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::particleTransferModels::dispersedBrownianDeposition::K(volScalarField& N)
const
{
    label i = Foam::flotationSystem::sectionNum(N);
    const scalar pi(constant::mathematical::pi);

    const volScalarField d(interface_.dispersed().d());
    const volScalarField D(model_[i].D());
    const volScalarField nu(interface_.continuous().thermo().nu());

    // Reynolds of the flow in the continuous phase induced by the presence of
    // the dispersed phase

    const volScalarField Re(interface_.Re());

    const volScalarField Sc(nu/D);
    const volScalarField Sh(Sh_->Sh(Re,Sc));

    return pi*Sh*d*D;
}

// ************************************************************************* //
