#include "LogMoMBreakupModel.H"
#include "phaseCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LogMoMBreakupModel, 0);
    defineRunTimeSelectionTable(LogMoMBreakupModel, dictionary);
}

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::LogMoMBreakupModel> Foam::LogMoMBreakupModel::New
(
    const diameterModels::LogMoM& logmom,
    const dictionary& dict
)
{
    word breakupModelType(dict.lookup("type"));

    Info<< "Selecting breakup model for " << breakupModelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(breakupModelType);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown breakup model "
            << breakupModelType << endl << endl
            << "Valid breakup models are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return ctorPtr(logmom, dict);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LogMoMBreakupModel::LogMoMBreakupModel
(
    const diameterModels::LogMoM& logmom,
    const dictionary& dict
)
:
    logmom_(logmom)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LogMoMBreakupModel::~LogMoMBreakupModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::phaseCompressibleTurbulenceModel&
Foam::LogMoMBreakupModel::continuousTurbulence() const
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
