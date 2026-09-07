#include "continuousTurbulentEntrainmentRate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace particleTransferModels
{
    defineTypeNameAndDebug(continuousTurbulentEntrainmentRate, 0);
    defineRunTimeSelectionTable(continuousTurbulentEntrainmentRate, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particleTransferModels::continuousTurbulentEntrainmentRate::continuousTurbulentEntrainmentRate
(
    continuousTurbulentEntrainment& transfer,
    const dictionary& dict
)
:
    transfer_(transfer)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::particleTransferModels::continuousTurbulentEntrainmentRate::
~continuousTurbulentEntrainmentRate()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::particleTransferModels::continuousTurbulentEntrainmentRate>
Foam::particleTransferModels::continuousTurbulentEntrainmentRate::New
(
    continuousTurbulentEntrainment& transfer,
    const dictionary& dict
)
{
    word modelType(dict.lookup("type"));

    Info<< "Selecting continuous entrainment rate" << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown model type " << modelType << endl << endl
            << "Valid model types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return ctorPtr(transfer, dict);
}

// ************************************************************************* //
