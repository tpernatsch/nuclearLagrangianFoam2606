#include "flotationSystem.H"
#include "fvm.H"
#include "fvc.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(flotationSystem, 0);
}

// * * * * * * * * * * * * * * * Local Functions  * * * * * * * * * * * * * //

namespace
{
    // Split an interface name of the form "<phase1>_<phase2>" (e.g.
    // "water_air") into references to its two constituent phaseModels,
    // by checking every ordered pair of the phaseSystem's actual phases -
    // this is the ESI-side replacement for Foundation's phaseInterface
    // word-parsing constructor, which ESI's phasePair doesn't have (it
    // always takes two phaseModel references directly).
    void splitInterfaceName
    (
        const Foam::phaseSystem& fluid,
        const Foam::word& interfaceName,
        const Foam::phaseModel*& phase1,
        const Foam::phaseModel*& phase2
    )
    {
        forAll(fluid.phases(), i)
        {
            forAll(fluid.phases(), j)
            {
                if (i == j)
                {
                    continue;
                }

                const Foam::word& name1 = fluid.phases()[i].name();
                const Foam::word& name2 = fluid.phases()[j].name();

                if (name1 + "_" + name2 == interfaceName)
                {
                    phase1 = &fluid.phases()[i];
                    phase2 = &fluid.phases()[j];
                    return;
                }
            }
        }

        FatalErrorInFunction
            << "Interface name " << interfaceName << " could not be split "
            << "into two known phases of the form <phase1>_<phase2>."
            << Foam::exit(Foam::FatalError);
    }

    Foam::phasePair makeInterfacePair
    (
        const Foam::phaseSystem& fluid,
        const Foam::word& interfaceName
    )
    {
        const Foam::phaseModel* phase1 = nullptr;
        const Foam::phaseModel* phase2 = nullptr;
        splitInterfaceName(fluid, interfaceName, phase1, phase2);
        return Foam::phasePair(*phase1, *phase2);
    }
}

// * * * * * * * * * * * * * * Static Member Functions* * * * * * * * * * * *//

Foam::label Foam::flotationSystem::sectionNum(const volScalarField& f)
{
    return atoi(f.name().substr(f.name().find(":") + 1).c_str()) - 1;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flotationSystem::flotationSystem
(
    const phaseSystem& fluid,
    const fvMesh& mesh,
    const word dictName
)
:
    IOdictionary
    (
        IOobject
        (
            dictName == "" ? "flotationProperties" : dictName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    fluid_(fluid),
    mesh_(mesh),
    rho_("rho", dimDensity, *this),
    nCorr_(lookupOrDefault<label>("nCorr", 2)),
    interface_(makeInterfacePair(fluid, word(lookup("interface")))),
    distribution_
    (
        particleSizeDistribution::New
        (
            subDict("distribution"),
            mesh
        ).ptr()
    ),
    particleModels_(2)
{
    // Set the particle model for each side of the interface

    particleModels_.set(0, new particleModel(*this, interface_.phase1()));
    particleModels_.set(1, new particleModel(*this, interface_.phase2()));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::flotationSystem::~flotationSystem()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::flotationSystem::solve()
{
    // Solve

    scalarList r0(distribution_->size(), 0.0);
    scalarList r(distribution_->size(), 0.0);
    labelList iters(distribution_->size(), 0);

    label iter;
    for (iter = 0; iter < nCorr_; iter++)
    {
        forAll(particleModels_, i)
        {
            List<solverPerformance> perfs(particleModels_[i].solve());

            if (iter == 0)
            {
                forAll(perfs, i)
                {
                    r0[i] = max(r0[i], perfs[i].initialResidual());
                }
            }

            if (iter == nCorr_-1)
            {
                forAll(perfs, i)
                {
                    r[i] = max(r[i], perfs[i].finalResidual());
                }
            }

            forAll(perfs, i)
            {
                iters[i] += perfs[i].nIterations();
            }
        }
    }

    Info<< typeName << ":" << endl;

    forAll(distribution_(), i)
    {
        Info<< "    Solving for section " << i+1
            << ", Initial residual = " << r0[i]
            << ", Final residual = " << r[i]
            << ", No Iterations = " << iters[i]/iter/2
            << ", No Corrections = " << iter << endl;
    }
}

bool Foam::flotationSystem::read()
{
    if (regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}

const Foam::phaseCompressibleTurbulenceModel&
Foam::flotationSystem::firstPhaseTurbulence() const
{
    return
        mesh_.lookupObject<phaseCompressibleTurbulenceModel>
        (
            IOobject::groupName
            (
                turbulenceModel::propertiesName,
                interface().phase1().name()
            )
        );
}

const Foam::phaseCompressibleTurbulenceModel&
Foam::flotationSystem::secondPhaseTurbulence() const
{
    return
        mesh_.lookupObject<phaseCompressibleTurbulenceModel>
        (
            IOobject::groupName
            (
                turbulenceModel::propertiesName,
                interface().phase2().name()
            )
        );
}

// ************************************************************************* //
