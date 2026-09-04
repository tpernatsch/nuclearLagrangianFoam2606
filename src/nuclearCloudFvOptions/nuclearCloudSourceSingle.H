/*---------------------------------------------------------------------------*\
Class
    Foam::fv::nuclearCloudSourceSingle

Description
    One-way coupled Lagrangian nuclear parcel cloud, run as an fvOption on
    top of a single-phase buoyant solver (buoyantBoussinesqPimpleFoam,
    buoyantBoussinesqSimpleFoam, ...) whose equations call
    fvOptions.correct(U) and provide real "U", "T" and "g" fields - i.e.
    the sibling of nuclearCloudSource that doesn't need a phaseSystem.

    Unlike nuclearCloudSourceSingle's carrier density/viscosity (which
    have no natural source in a Boussinesq solver - it only ever solves
    for kinematic quantities), "T" and "g" are both genuine persistent
    objects the solver already owns, so - unlike phaseModel::rho()/mu()
    in nuclearCloudSource - there is no tmp<> dangling-reference trap to
    work around for them. T is looked up as a real MUST_READ field on the
    mesh registry. g is NOT: in this OpenFOAM version
    readGravitationalAcceleration.H gets it via
    meshObjects::gravity::New(runTime), a MeshObject cached on the Time
    database rather than a plain field registered on the mesh - a
    mesh_.lookupObject<uniformDimensionedVectorField>("g") fails ("failed
    lookup of g ... available objects of type
    uniformDimensionedVectorField: 0()") even though the solver clearly
    read it at startup. Fetch it the same way the solver does instead:
    meshObjects::gravity::New(mesh_.time()) (gravityMeshObject.H) -
    Foam::meshObjects::gravity publicly derives from
    uniformDimensionedVectorField, so it's usable anywhere a
    dimensionedVector is expected. Only rho/mu still need owned copies,
    for the same reason as kinematicCloudSource: the cloud requires real
    volScalarFields, not dimensionedScalar constants.

    Never touches the momentum/temperature matrices: addSup() is left at
    the fv::option default (a no-op). One-way coupling only.

    Example usage, in constant/fvOptions:
    \verbatim
    FPCloud
    {
        type                nuclearCloudSourceSingle;

        nuclearCloudSourceSingleCoeffs
        {
            cloudName       FPCloud;  // optional, defaults to <name>Cloud
            rho             1.2;      // carrier density [kg/m3]
            nu              1.5e-05;  // carrier kinematic viscosity [m2/s]
        }
    }
    \endverbatim

SourceFiles
    nuclearCloudSourceSingle.C

\*---------------------------------------------------------------------------*/

#ifndef nuclearCloudSourceSingle_H
#define nuclearCloudSourceSingle_H

#include "fvOption.H"
#include "basicNuclearCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class nuclearCloudSourceSingle Declaration
\*---------------------------------------------------------------------------*/

class nuclearCloudSourceSingle
:
    public option
{
    // Private Data

        //- Name of the cloud
        word cloudName_;

        //- Uniform carrier density, held constant
        dimensionedScalar rho0_;

        //- Uniform carrier dynamic viscosity (= rho0_ * nu)
        dimensionedScalar mu0_;

        //- Persistent carrier density field (see class doc - required by
        //  the cloud constructor's type, not to dodge a dangling tmp<>)
        volScalarField rhoc_;

        //- Persistent carrier dynamic viscosity field
        volScalarField muc_;

        //- Persistent vorticity of the carrier velocity, refreshed every
        //  correct() call
        volVectorField vort_;

        //- The cloud itself; created lazily on the first correct() call
        autoPtr<basicNuclearCloud> cloud_;


public:

    //- Runtime type information
    TypeName("nuclearCloudSourceSingle");


    // Constructors

        nuclearCloudSourceSingle
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        nuclearCloudSourceSingle(const nuclearCloudSourceSingle&) = delete;

        //- No copy assignment
        void operator=(const nuclearCloudSourceSingle&) = delete;


    //- Destructor
    virtual ~nuclearCloudSourceSingle()
    {}


    // Member Functions

        // Evaluation

            //- Evolve the cloud. Triggered once per PIMPLE iteration by
            //  fvOptions.correct(U). Deliberately does not touch the
            //  fvMatrix - one-way coupling only.
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
