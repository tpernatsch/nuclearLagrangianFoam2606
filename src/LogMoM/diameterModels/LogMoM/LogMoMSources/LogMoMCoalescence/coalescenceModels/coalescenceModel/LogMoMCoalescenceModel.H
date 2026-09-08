#ifndef LogMoMCoalescenceModel_H
#define LogMoMCoalescenceModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "LogMoM.H"
#include "runTimeSelectionTables.H"
#include "phaseCompressibleTurbulenceModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class LogMoMCoalescenceModel Declaration
\*---------------------------------------------------------------------------*/

class LogMoMCoalescenceModel
{
protected:

    // Protected data

        //- Reference to the LogMoM diameter model
        const diameterModels::LogMoM& logmom_;

public:

    //- Runtime type information
    TypeName("coalescenceModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            LogMoMCoalescenceModel,
            dictionary,
            (
                const diameterModels::LogMoM& logmom,
                const dictionary& dict
            ),
            (logmom, dict)
        );


    // Constructors

        // Construct from ditionary
        LogMoMCoalescenceModel
        (
            const diameterModels::LogMoM& logmom,
            const dictionary& dict
        );


    //- Destructor
    virtual ~LogMoMCoalescenceModel();


    // Selectors

        static autoPtr<LogMoMCoalescenceModel> New
        (
            const diameterModels::LogMoM& logmom,
            const dictionary& dict
        );


    // Member Functions

        // Return a field of coalescence rates given two particle sizes
        virtual tmp<volScalarField> rate
        (
            const volScalarField& di,
            const volScalarField& dj
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
