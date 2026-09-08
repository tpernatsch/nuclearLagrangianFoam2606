#include "GaussQuadrature.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GaussQuadrature, 0);
    defineRunTimeSelectionTable(GaussQuadrature, dictionary);
}

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::GaussQuadrature> Foam::GaussQuadrature::New
(
    const word quadrature,
    const label numNodes
)
{
    const word GaussQuadratureType(quadrature+Foam::name(numNodes));

    auto* ctorPtr = dictionaryConstructorTable(GaussQuadratureType);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown GaussQuadratureType type "
            << GaussQuadratureType << endl << endl
            << "Valid GaussQuadrature types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return ctorPtr();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GaussQuadrature::GaussQuadrature()
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GaussQuadrature::~GaussQuadrature()
{}

// ************************************************************************* //
