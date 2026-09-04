/*---------------------------------------------------------------------------*\
Class
    Foam::fv::nuclearCloudSource

Description
    One-way coupled Lagrangian nuclear parcel cloud, run as an fvOption
    on top of a phase of reactingMultiphaseEulerFoam (or any solver whose
    equations call fvOptions.correct(U.<phase>)).

    The cloud reads U, T, rho and mu from the named carrier phase and
    evolves (advection, decay, decay heat bookkeeping) once per outer
    (PIMPLE) iteration. It never adds anything to the phase momentum or
    energy matrices: addSup() is left at the fv::option default (a no-op),
    so no new solver is needed and no two-way coupling terms exist to
    remove - there simply are none.

    Example usage, in constant/fvOptions:
    \verbatim
    FPCloud
    {
        type                nuclearCloudSource;

        nuclearCloudSourceCoeffs
        {
            phase           salt;
            cloudName       FPCloud;    // optional, defaults to <name>Cloud
        }
    }
    \endverbatim

    The named phase must exist in constant/phaseProperties, and the case
    must provide FPCloudProperties (particle injection, forces, etc.) the
    same way a nuclearParcelFoam case does.

SourceFiles
    nuclearCloudSource.C

\*---------------------------------------------------------------------------*/

#ifndef nuclearCloudSource_H
#define nuclearCloudSource_H

#include "fvOption.H"
#include "phaseSystem.H"
#include "phaseModel.H"
#include "basicNuclearCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class nuclearCloudSource Declaration
\*---------------------------------------------------------------------------*/

class nuclearCloudSource
:
    public option
{
    // Private Data

        //- Name of the carrier phase this cloud rides on
        word phaseName_;

        //- Name of the cloud
        word cloudName_;

        //- Persistent copy of the carrier density.
        //  phaseModel::rho() returns a tmp<>, but NuclearCloud stores its
        //  carrier fields by reference for the cloud's whole lifetime, so
        //  a tmp<> cannot be handed to it directly - the referenced field
        //  would be destroyed as soon as the tmp<> is, leaving a dangling
        //  reference. rhoc_/muc_/vort_ are real, persistent fields owned
        //  by this option and refreshed every call, so the cloud always
        //  has something valid to point at.
        volScalarField rhoc_;

        //- Persistent copy of the carrier dynamic viscosity (same reason)
        volScalarField muc_;

        //- Persistent vorticity of the carrier velocity, needed by the
        //  nuclear parcel's thermophoresis force
        volVectorField vort_;

        //- The cloud itself; created lazily on the first correct() call,
        //  once all phase fields are guaranteed to exist
        autoPtr<basicNuclearCloud> cloud_;


    // Private Member Functions

        //- Return the carrier phase looked up from the mesh's phaseSystem
        phaseModel& phase() const;

        //- Refresh rhoc_/muc_/vort_ from the current phase fields, and
        //  build the cloud on the first call
        void updateCarrierFields();


public:

    //- Runtime type information
    TypeName("nuclearCloudSource");


    // Constructors

        nuclearCloudSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        nuclearCloudSource(const nuclearCloudSource&) = delete;

        //- No copy assignment
        void operator=(const nuclearCloudSource&) = delete;


    //- Destructor
    virtual ~nuclearCloudSource()
    {}


    // Member Functions

        // Evaluation

            //- Evolve the cloud. Triggered once per phase, per PIMPLE
            //  iteration, by fvOptions.correct(U.<phase>). Deliberately
            //  does not touch the fvMatrix - one-way coupling only.
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
