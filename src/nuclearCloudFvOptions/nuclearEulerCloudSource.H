/*---------------------------------------------------------------------------*\
Class
    Foam::fv::nuclearEulerCloudSource

Description
    One-way coupled, two-phase-Euler-aware Lagrangian nuclear parcel cloud,
    run as an fvOption on top of two phases of reactingMultiphaseEulerFoam
    (or any Euler-Euler solver whose equations call
    fvOptions.correct(U.<phase>)).

    This is the two-phase counterpart of nuclearCloudSource: rather than
    feeding a plain basicNuclearCloud from a single named phase, it feeds a
    basicNuclearEulerCloud from *two* named phases -- the "primary" phase
    the parcels are actually tracked/advected through (drag, buoyancy,
    patch interaction all act against this phase, unchanged from
    nuclearCloudSource), and a "secondary" phase whose volume fraction,
    velocity, density, viscosity, temperature and turbulence closure
    fields are additionally exposed to the cloud for the two-phase heat
    transfer and turbulent dispersion closures described in
    NuclearEulerCloud.H / NuclearEulerParcel.H.

    Exactly like nuclearCloudSource, it never adds anything to either
    phase's momentum or energy matrices: addSup() is left at the
    fv::option default (a no-op), so no new solver is needed.

    Example usage, in constant/fvOptions:
    \verbatim
    FPCloud
    {
        type                nuclearEulerCloudSource;

        nuclearEulerCloudSourceCoeffs
        {
            phase           liquid;     // primary/tracking phase
            phase2          gas;        // secondary phase
            cloudName       FPCloud;    // optional, defaults to <name>Cloud
        }
    }
    \endverbatim

    Both named phases must exist in constant/phaseProperties, and the case
    must provide FPCloudProperties (particle injection, forces, etc.) the
    same way a nuclearParcelFoam case does.

SourceFiles
    nuclearEulerCloudSource.C

\*---------------------------------------------------------------------------*/

#ifndef nuclearEulerCloudSource_H
#define nuclearEulerCloudSource_H

#include "fvOption.H"
#include "phaseSystem.H"
#include "phaseModel.H"
#include "basicNuclearEulerCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                  Class nuclearEulerCloudSource Declaration
\*---------------------------------------------------------------------------*/

class nuclearEulerCloudSource
:
    public option
{
    // Private Data

        //- Name of the primary (tracking) carrier phase
        word phaseName_;

        //- Name of the secondary Euler carrier phase
        word phase2Name_;

        //- Name of the cloud
        word cloudName_;

        // Persistent copies of fields phaseModel returns as tmp<>

            //- Primary phase volume fraction
            volScalarField alpha1c_;

            //- Primary phase density
            volScalarField rhoc_;

            //- Primary phase dynamic viscosity
            volScalarField muc_;

            //- Primary phase turbulence kinetic energy
            volScalarField k1c_;

            //- Primary phase turbulence dissipation rate
            volScalarField epsilon1c_;

            //- Secondary phase volume fraction
            volScalarField alpha2c_;

            //- Secondary phase density
            volScalarField rho2c_;

            //- Secondary phase dynamic viscosity
            volScalarField mu2c_;

            //- Secondary phase turbulence kinetic energy
            volScalarField k2c_;

            //- Secondary phase turbulence dissipation rate
            volScalarField epsilon2c_;

        //- The cloud itself; created lazily on the first correct() call,
        //  once all phase fields are guaranteed to exist
        autoPtr<basicNuclearEulerCloud> cloud_;


    // Private Member Functions

        //- Return the named phase looked up from the mesh's phaseSystem
        phaseModel& phaseByName(const word& phaseName) const;

        //- Return the primary (tracking) phase
        phaseModel& phase() const;

        //- Return the secondary phase
        phaseModel& phase2() const;

        //- Return the named phase's turbulence model, looked up from the
        //  mesh object registry
        const phaseCompressibleTurbulenceModel& turbulenceByName
        (
            const word& phaseName
        ) const;

        //- Return true if the named phase's turbulence model is RAS
        bool isRASByName(const word& phaseName) const;

        //- Refresh every persistent copy from the current phase fields,
        //  and build the cloud on the first call
        void updateCarrierFields();


public:

    //- Runtime type information
    TypeName("nuclearEulerCloudSource");


    // Constructors

        nuclearEulerCloudSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        nuclearEulerCloudSource(const nuclearEulerCloudSource&) = delete;

        //- No copy assignment
        void operator=(const nuclearEulerCloudSource&) = delete;


    //- Destructor
    virtual ~nuclearEulerCloudSource()
    {}


    // Member Functions

        // Evaluation

            //- Evolve the cloud, one-way coupling only
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
