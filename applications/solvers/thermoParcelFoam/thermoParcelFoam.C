#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "surfaceFilmModel.H"
#include "basicKinematicCloud.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "syncTools.H"

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow"
        " with kinematic particle clouds"
        " and surface film modelling. For the thermic part the"
        " Boussinesq approximation is used."
    );

    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createRegionControls.H"
    #include "createUfIfPresent.H"

    turbulence->validate();

    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setMultiRegionDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        FP1.storeGlobalPositions();
        FP1.evolve();

        Foam::label totalCells = mesh.nCells();
        for (Foam::label cellI = 0; cellI < totalCells; ++cellI)
        {
            const Foam::DynamicList<Foam::KinematicParcel<Foam::particle>*>& particlesInCell = FP1.cellOccupancy()[cellI];
            scalar nParticles = 0;
            forAll(particlesInCell, parcelI)
            {
                const Foam::KinematicParcel<Foam::particle>& FP1 = *particlesInCell[parcelI];
                nParticles += FP1.nParticle();
            }
            scalar cellVolume = mesh.V()[cellI];
            Co[cellI] = (cellVolume > VSMALL) ? nParticles / cellVolume : 0.0;
        }


        Co.setOriented();
        Co.correctBoundaryConditions();
        laplacianCo = fvc::laplacian(dimensionedScalar("one", dimless, 1.0), Co);

        surfaceFilm.evolve();

        if (solvePrimaryRegion)
        {
            while (pimple.loop())
            {
                #include "UEqn.H"
                #include "TEqn.H"

                while (pimple.correct())
                {
                    #include "pEqn.H"
                }

                if (pimple.turbCorr())
                {
                    laminarTransport.correct();
                    turbulence->correct();
                }
            }
        }

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;
    return 0;
}
