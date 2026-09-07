#include "particleTransferModel.H"
#include "particleModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(particleTransferModel, 0);
    defineRunTimeSelectionTable(particleTransferModel, dictionary);
}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::particleTransferModel> Foam::particleTransferModel::New
(
    const word& type,
    const particleModel& model,
    const dictionary& dict
)
{
    auto* ctorPtr = dictionaryConstructorTable(type);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown particle transfer model type "
            << type << nl << nl
            << "Valid particle transfer model types : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<particleTransferModel>(ctorPtr(model, dict));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
