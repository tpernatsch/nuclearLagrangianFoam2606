#ifndef LogMoMBreakupModel_H
#define LogMoMBreakupModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LogMoM.H"
#include "runTimeSelectionTables.H"
#include "phaseCompressibleTurbulenceModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class LogMoMBreakupModel Declaration
\*---------------------------------------------------------------------------*/

class LogMoMBreakupModel
{
protected:

    // Protected data

        //- Reference to the LogMoM diameter model
        const diameterModels::LogMoM& logmom_;

public:

    //- Runtime type information
    TypeName("breakupModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            LogMoMBreakupModel,
            dictionary,
            (
                const diameterModels::LogMoM& logmom,
                const dictionary& dict
            ),
            (logmom, dict)
        );


    // Constructors

        // Construct from ditionary
        LogMoMBreakupModel
        (
            const diameterModels::LogMoM& logmom,
            const dictionary& dict
        );


    //- Destructor
    virtual ~LogMoMBreakupModel();


    // Selectors

        static autoPtr<LogMoMBreakupModel> New
        (
            const diameterModels::LogMoM& logmom,
            const dictionary& dict
        );


    // Member Functions

        // Return a field of binary breakup rates given two particle diameters
        virtual tmp<volScalarField> binaryRate
        (
            const volScalarField& d1,
            const volScalarField& d2
        ) const = 0;

        // Return the continuous phase turbulence model
        virtual const phaseCompressibleTurbulenceModel&
            continuousTurbulence() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
