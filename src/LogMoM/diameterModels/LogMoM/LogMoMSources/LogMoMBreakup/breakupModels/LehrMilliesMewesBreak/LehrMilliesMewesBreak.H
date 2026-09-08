#ifndef LehrMilliesMewes_H
#define LehrMilliesMewes_H

#include "LogMoMBreakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace breakupModels
{

/*---------------------------------------------------------------------------*\
                    Class LehrMilliesMewes Declaration
\*---------------------------------------------------------------------------*/

class LehrMilliesMewes
:
    public LogMoMBreakupModel
{
public:

    //- Runtime type information
    TypeName("LehrMilliesMewes");

    // Constructor

        LehrMilliesMewes
        (
            const diameterModels::LogMoM& logmom,
            const dictionary& dict
        );


    //- Destructor
    virtual ~LehrMilliesMewes()
    {}


    // Member Functions

        // Return a field of binary breakup rates given two particle diameters
        virtual tmp<volScalarField> binaryRate
        (
            const volScalarField& d1,
            const volScalarField& d2
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace breakupModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
