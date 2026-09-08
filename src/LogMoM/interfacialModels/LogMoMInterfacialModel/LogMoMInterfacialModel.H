#ifndef LogMoMInterfacialModel_H
#define LogMoMInterfacialModel_H

#include "LogMoM.H"
#include "GaussQuadrature.H"
#include "phasePair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class LogMoMInterfacialModel;

/*---------------------------------------------------------------------------*\
                    Class LogMoMInterfacialModel Declaration
\*---------------------------------------------------------------------------*/

class LogMoMInterfacialModel
{
public:

    // Public data types

        //- Integration mode enumeration
        enum class integrationMode
        {
            none,
            quadrature,
            diameter
        };

        //- Integration mode names
        static const Enum<integrationMode> integrationModeNames_;

private:

    // Private Data

        //- Integration mode
        integrationMode integrationMode_;

        //- For the diameter integration mode we need the p and q coefficients
        //  of the average diameter d_{p,q}, which are provided as polynomials
        //  of gamma with the first entry multiplying gamma^0, the second entry
        //  gamma^1, etc.
        scalarList p_;
        scalarList q_;

        //- By default we have p = {0,1} and q = {2,0}, i.e., d_{gamma,2}, so
        //  that the Sauter-mean diameter is obtained for gamma = 3.
        static const scalarList pd_;
        static const scalarList qd_;

        //- The dispersed-phase-ordered interface. Bound to the same phasePair
        //  object handed to the owning dragModel/liftModel, which the
        //  phaseSystem constructs already ordered (dispersed, continuous)
        //  for a "x_dispersedIn_y" entry, so dispersed()/continuous() resolve
        //  correctly through phasePair's own virtual dispatch.
        const phasePair& pair_;

        //- Gauss-Hermite quadrature rules
        autoPtr<GaussQuadrature> quadrature_;


public:

    //- Runtime type information
    TypeName("LogMoM");


    // Constructors

        //- Construct from a dictionary and an interface
        LogMoMInterfacialModel
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~LogMoMInterfacialModel();


    // Member Functions

        //- Return the LogMoM model corresponding to the interface
        const diameterModels::LogMoM& logmom() const;

        //- Return the current moment order to be integrated. This is retrieved
        //  from the LogMoM model.
        scalar gamma() const;

        //- Evaluate the function (1/M_gamma) int [d^beta*f()*n(d)] dd on
        //  [0,inf], where f() is provided as a function pointer. If the
        //  integration mode is 'quadrature', centered Gauss-Hermite quadrature
        //  is used. If the integration mode is 'diameter', a specified mean
        //  diameter d_{p(beta),q(beta)} is used as an integration point to
        //  evaluate the integral.
        template
        <
            class Type,
            template<class> class PatchField,
            class GeoMesh,
            class Model
        >
        tmp<GeometricField<Type, PatchField, GeoMesh>> evaluate
        (
            const scalar beta,
            const scalar gamma,
            tmp<GeometricField<Type, PatchField, GeoMesh>> (Model::*f)() const,
            const Model& model
        ) const;

        // Helper function for when beta equals gamma
        template
        <
            class Type,
            template<class> class PatchField,
            class GeoMesh,
            class Model
        >
        tmp<GeometricField<Type, PatchField, GeoMesh>> evaluate
        (
            const scalar gamma,
            tmp<GeometricField<Type, PatchField, GeoMesh>> (Model::*f)() const,
            const Model& model
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LogMoMInterfacialModelTemplates.C"

#endif

// ************************************************************************* //
