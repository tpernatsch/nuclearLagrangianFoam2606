#include "LogMoMDragModel.H"
#include "addToRunTimeSelectionTable.H"
#include "swarmCorrection.H"
#include "noSwarm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(LogMoMDragModel, 0);
    addToRunTimeSelectionTable(dragModel, LogMoMDragModel, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::LogMoMDragModel::LogMoMDragModel
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    dragModel(pair, registerObject),
    LogMoMInterfacialModel(dict, pair),
    dragModelPtr_
    (
        dragModel::dictionaryConstructorTable
        (
            dict.subDict("drag").get<word>("type")
        )
        (
            dict.subDict("drag"),
            pair,
            false
        )
    ),
    swarmCorrection_
    (
        dict.found("swarmCorrection")
      ? swarmCorrection::New(dict.subDict("swarmCorrection"), pair).ptr()
      : new swarmCorrections::noSwarm(dict, pair)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::LogMoMDragModel::~LogMoMDragModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::LogMoMDragModel::CdRe() const
{
    return dragModelPtr_->CdRe();
}

Foam::tmp<Foam::volScalarField> Foam::dragModels::LogMoMDragModel::Ki() const
{
    // We need to take the gamma'th moment of Cd*Re/d^2. So, instead, we take
    // the (gamma - 2)'th moment of Cd*Re, so that the centered Gauss-Hermite
    // integration is more accurate.

    return
        0.75
      * this->evaluate(gamma() - 2, gamma(), &LogMoMDragModel::CdRe, *this)
      * swarmCorrection_->Cs()
      * dragModel::pair_.continuous().rho()
      * dragModel::pair_.continuous().nu();
}

// ************************************************************************* //
