#include "continuousParticleTransferModel.H"
#include "particleModel.H"
#include "flotationSystem.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particleTransferModels::continuousParticleTransferModel::
continuousParticleTransferModel
(
    const particleModel& model,
    const dictionary& dict
)
:
    particleTransferModel(model),
    interface_(model.phase(), model.otherPhase())
{
    const phaseSystem& fluid = model.system().fluid();

    const dictionary& blendingDict =
        fluid.subDict("blending").found("particleTransfer")
      ? fluid.subDict("blending").subDict("particleTransfer")
      : fluid.subDict("blending").subDict("default");

    blending_.set
    (
        blendingMethod::New
        (
            "particleTransfer",
            blendingDict,
            wordList
            {
                interface_.phase1().name(),
                interface_.phase2().name()
            }
        ).ptr()
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::particleTransferModels::continuousParticleTransferModel::R
(
    volScalarField& N
) const
{
    // Only activate when the first phase (i.e. the phase to which the particle
    // model belongs) is considered as dispersed

    const volScalarField blending
    (
        blending_->f1(interface_.dispersed(), interface_.continuous())
    );

    return -fvm::Sp(blending*this->K(N), N);
}

// ************************************************************************* //
