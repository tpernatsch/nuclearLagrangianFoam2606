#include "LogMoMCoalescenceModel.H"
#include "phaseCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LogMoMCoalescenceModel, 0);
    defineRunTimeSelectionTable(LogMoMCoalescenceModel, dictionary);
}

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::LogMoMCoalescenceModel> Foam::LogMoMCoalescenceModel::New
(
    const diameterModels::LogMoM& logmom,
    const dictionary& dict
)
{
    word coalescenceModelType(dict.lookup("type"));

    Info<< "Selecting coalescence model for " << coalescenceModelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(coalescenceModelType);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown coalescence model "
            << coalescenceModelType << endl << endl
            << "Valid coalescence models are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return ctorPtr(logmom, dict);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LogMoMCoalescenceModel::LogMoMCoalescenceModel
(
    const diameterModels::LogMoM& logmom,
    const dictionary& dict
)
:
    logmom_(logmom)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LogMoMCoalescenceModel::~LogMoMCoalescenceModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::phaseCompressibleTurbulenceModel&
Foam::LogMoMCoalescenceModel::continuousTurbulence() const
{
    return
        logmom_.phase().mesh().lookupObject<phaseCompressibleTurbulenceModel>
        (
            IOobject::groupName
            (
                turbulenceModel::propertiesName,
                logmom_.continuousPhase().name()
            )
        );
}

// ************************************************************************* //
