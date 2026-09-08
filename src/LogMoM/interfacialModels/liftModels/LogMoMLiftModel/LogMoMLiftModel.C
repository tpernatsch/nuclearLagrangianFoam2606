#include "LogMoMLiftModel.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcCurl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace liftModels
{
    defineTypeNameAndDebug(LogMoMLiftModel, 0);
    addToRunTimeSelectionTable(liftModel, LogMoMLiftModel, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liftModels::LogMoMLiftModel::LogMoMLiftModel
(
    const dictionary& dict,
    const phasePair& pair
)
:
    liftModel(dict, pair),
    LogMoMInterfacialModel(dict, pair),
    liftModelPtr_
    (
        liftModel::New(dict.subDict("lift"), pair).ptr()
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::liftModels::LogMoMLiftModel::~LogMoMLiftModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::liftModels::LogMoMLiftModel::Cl() const
{
    return liftModelPtr_->Cl();
}

Foam::tmp<Foam::volVectorField> Foam::liftModels::LogMoMLiftModel::Fi() const
{
    return
        this->evaluate(gamma(), &LogMoMLiftModel::Cl, *this)
      * liftModel::pair_.continuous().rho()
      * (
            liftModel::pair_.Ur() ^ fvc::curl(liftModel::pair_.continuous().U())
        );
}

// ************************************************************************* //
