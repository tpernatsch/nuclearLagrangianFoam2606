/*---------------------------------------------------------------------------*\
Class
    Foam::fv::kinematicCloudSource

Description
    One-way coupled Lagrangian kinematic parcel cloud, run as an fvOption
    on top of any single-phase incompressible solver whose equations call
    fvOptions.correct(U) (pimpleFoam, simpleFoam, icoFoam, ...).

    This is the single-phase counterpart of nuclearCloudSource: instead of
    looking a carrier phase up from a phaseSystem, it reads a constant
    carrier density/viscosity from its own dictionary entries and rides
    directly on the solver's own "U" field. It never adds anything to the
    momentum equation: addSup() is left at the fv::option default (a
    no-op), so there is nothing to configure for one-way coupling - it
    simply never touches the matrix.

    Example usage, in constant/fvOptions:
    \verbatim
    particleCloud
    {
        type                kinematicCloudSource;

        kinematicCloudSourceCoeffs
        {
            cloudName       kinematicCloud; // optional, defaults to <name>Cloud
            rho             1000;           // carrier density [kg/m3]
            nu              1e-05;          // carrier kinematic viscosity [m2/s]
        }
    }
    \endverbatim

    The case must also provide a <cloudName>Properties dictionary (particle
    injection, forces, patch interaction, etc.), the same way any
    Lagrangian solver's case does.

SourceFiles
    kinematicCloudSource.C

\*---------------------------------------------------------------------------*/

#ifndef kinematicCloudSource_H
#define kinematicCloudSource_H

#include "fvOption.H"
#include "basicKinematicCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class kinematicCloudSource Declaration
\*---------------------------------------------------------------------------*/

class kinematicCloudSource
:
    public option
{
    // Private Data

        //- Name of the cloud
        word cloudName_;

        //- Uniform carrier density, held constant for this simple carrier
        dimensionedScalar rho0_;

        //- Uniform carrier dynamic viscosity (= rho0_ * nu), held constant
        dimensionedScalar mu0_;

        //- Gravity used by the cloud's particle forces (defaults to zero
        //  if the case has no constant/g and none is given here)
        dimensionedVector g_;

        //- Persistent carrier density field.
        //  KinematicCloud stores its carrier fields by reference for the
        //  cloud's whole lifetime (see nuclearCloudSource's note on the
        //  same trap with phaseModel::rho()). rho0_/mu0_ are plain
        //  dimensionedScalar constants here, not tmp<> fields, so there
        //  is no dangling-reference risk from them directly - but the
        //  cloud constructor still requires a real volScalarField (it
        //  may vary spatially in general), so rhoc_/muc_ exist purely to
        //  satisfy that type, initialised once and never reassigned.
        volScalarField rhoc_;

        //- Persistent carrier dynamic viscosity field (see rhoc_)
        volScalarField muc_;

        //- Persistent vorticity of the carrier velocity, refreshed every
        //  correct() call since U changes every step
        volVectorField vort_;

        //- The cloud itself; created lazily on the first correct() call
        autoPtr<basicKinematicCloud> cloud_;


public:

    //- Runtime type information
    TypeName("kinematicCloudSource");


    // Constructors

        kinematicCloudSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        kinematicCloudSource(const kinematicCloudSource&) = delete;

        //- No copy assignment
        void operator=(const kinematicCloudSource&) = delete;


    //- Destructor
    virtual ~kinematicCloudSource()
    {}


    // Member Functions

        // Evaluation

            //- Evolve the cloud. Triggered once per PIMPLE/SIMPLE
            //  iteration by fvOptions.correct(U). Deliberately does not
            //  touch the fvMatrix - one-way coupling only.
            virtual void correct(volVectorField& U);


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
