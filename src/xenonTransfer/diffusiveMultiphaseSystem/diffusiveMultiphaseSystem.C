/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "diffusiveMultiphaseSystem.H"

#include "fixedValueFvsPatchFields.H"
#include "Time.H"
#include "subCycle.H"
#include "fvcMeshPhi.H"

#include "surfaceInterpolate.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcDiv.H"
#include "fvcDdt.H"
#include "fvcFlux.H"
#include "fvmDdt.H"
#include "fvcAverage.H"
#include "fvMatrix.H"
#include "fvmSup.H"
#include "CMULES.H"

namespace Foam
{
namespace multiphaseInter
{

defineTypeNameAndDebug(diffusiveMultiphaseSystem, 0);


diffusiveMultiphaseSystem::diffusiveMultiphaseSystem
(
    const fvMesh& mesh
)
:
    multiphaseSystem(mesh)
{}


void diffusiveMultiphaseSystem::solve()
{
    const dictionary& alphaControls = mesh_.solverDict("alpha");
    const label nAlphaSubCycles(alphaControls.get<label>("nAlphaSubCycles"));

    const scalar totalDeltaTValue = mesh_.time().deltaTValue();

    PtrList<volScalarField> alphaStart(phases_.size());

    label phasei = 0;
    for (const phaseModel& phase : phases_)
    {
        alphaStart.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "alphaStart." + phase.name(),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                phase
            )
        );

        ++phasei;
    }

    const Time& runTime = mesh_.time();

    if (nAlphaSubCycles > 1)
    {
        surfaceScalarField rhoPhiSum
        (
            IOobject
            (
                "rhoPhiSum",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", rhoPhi_.dimensions(), 0.0)
        );

        PtrList<surfaceScalarField> alphaPhiSums(phases_.size());

        label phasei = 0;
        for (phaseModel& phase : phases_)
        {
            alphaPhiSums.set
            (
                phasei,
                new surfaceScalarField
                (
                    IOobject
                    (
                        "alphaPhiSum." + phase.name(),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar
                    (
                        "zero",
                        phase.alphaPhi().dimensions(),
                        0.0
                    )
                )
            );

            ++phasei;
        }

        const scalar totalDeltaT = mesh_.time().deltaTValue();

        for
        (
            subCycleTime alphaSubCycle
            (
                const_cast<Time&>(runTime),
                nAlphaSubCycles
            );
            !(++alphaSubCycle).end();
        )
        {
            solveAlphas();

            const scalar w = mesh_.time().deltaTValue()/totalDeltaT;

            rhoPhiSum += w*rhoPhi_;

            phasei = 0;
            for (phaseModel& phase : phases_)
            {
                alphaPhiSums[phasei] += w*phase.alphaPhi();
                ++phasei;
            }
        }

        rhoPhi_ = rhoPhiSum;

        phasei = 0;
        for (phaseModel& phase : phases_)
        {
            phase.alphaPhi() = alphaPhiSums[phasei];
            ++phasei;
        }

        static label lastPrintedNAlphaSubCycles = -1;

        if (lastPrintedNAlphaSubCycles != nAlphaSubCycles)
        {
            Info<< "diffusiveMultiphaseSystem: using time-averaged phase.alphaPhi() over nAlphaSubCycles="
                << nAlphaSubCycles << nl;

            lastPrintedNAlphaSubCycles = nAlphaSubCycles;
        }
    }
    else
    {
        solveAlphas();
    }
}

} // End namespace multiphaseInter
} // End namespace Foam

// ************************************************************************* //