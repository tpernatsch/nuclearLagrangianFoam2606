#ifndef LogMoMLiftModel_H
#define LogMoMLiftModel_H

#include "liftModel.H"
#include "LogMoMInterfacialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace liftModels
{

/*---------------------------------------------------------------------------*\
                       Class LogMoMLiftModel Declaration
\*---------------------------------------------------------------------------*/

class LogMoMLiftModel
:
    public liftModel,
    public LogMoMInterfacialModel
{
    // Private Data

        //- The lift model
        autoPtr<liftModel> liftModelPtr_;

public:

    //- Runtime type information
    TypeName("LogMoM");


    // Constructors

        //- Construct from a dictionary and an interface
        LogMoMLiftModel
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~LogMoMLiftModel();


    // Member Functions

        //- Return lift coefficient
        virtual tmp<volScalarField> Cl() const;

        //- Return phase-intensive lift force
        virtual tmp<volVectorField> Fi() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace liftModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
