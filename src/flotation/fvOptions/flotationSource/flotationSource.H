/*---------------------------------------------------------------------------*\
Class
    Foam::fv::flotationSource

Description
    ESI v2606 port of flotationFoam's Foam::fv::flotation (an OpenFOAM-13
    fvModel). Runs the flotation particle-population transport equations
    (sectional N.<phase>:<section> fields) on top of a stock
    reactingMultiphaseEulerFoam case, exactly the way flotationFoam bolts
    onto Foundation's multiphaseEuler solver module: as a runtime-selected
    library, loaded via libs and selected from constant/fvOptions, with no
    new solver executable.

    Unlike the nuclear/kinematic cloud fvOptions in this same tree
    (nuclearCloudFvOptions, kinematicCloudFvOptions), this one does not
    override addSup() either - flotationSystem's own particleModel::solve()
    solves its N equations directly (see particleModel.C), the same way it
    always has. correct() is purely the per-timestep "run the solve" hook,
    mirroring fv::flotation::correct() calling system_->solve().

    Example usage, in constant/fvOptions:
    \verbatim
    flotation
    {
        type    flotationSource;
    }
    \endverbatim

    Parameters are read from constant/flotationProperties (or the file
    named by the optional "dict" entry), exactly as in the original
    flotationFoam.

SourceFiles
    flotationSource.C

\*---------------------------------------------------------------------------*/

#ifndef flotationSource_H
#define flotationSource_H

#include "fvOption.H"
#include "flotationSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class flotationSource Declaration
\*---------------------------------------------------------------------------*/

class flotationSource
:
    public option
{
    // Private Data

        //- The flotation system itself
        autoPtr<flotationSystem> system_;


public:

    //- Runtime type information
    TypeName("flotationSource");


    // Constructors

        flotationSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        flotationSource(const flotationSource&) = delete;

        //- No copy assignment
        void operator=(const flotationSource&) = delete;


    //- Destructor
    virtual ~flotationSource()
    {}


    // Member Functions

        // Evaluation

            //- Run the flotation system's own solve() once per timestep.
            //  Triggered off correct(U.<phase>) for the first (arbitrary)
            //  phase field this option is registered against, purely as a
            //  once-per-outer-iteration clock tick - addSup() is left at
            //  the fv::option default (a no-op): flotationSystem solves
            //  its own N equations directly, exactly as upstream.
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
