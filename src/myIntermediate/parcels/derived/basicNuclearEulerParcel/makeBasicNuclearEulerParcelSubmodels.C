/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tommaso Pernatsch
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

Description
    Instantiates the same submodel families as basicNuclearParcel. Every
    one of these submodel classes is templated purely on the CloudType and
    calls back into it only through the interface already provided by
    KinematicCloud/NuclearCloud (forces(), dispersion(), heatTransfer(),
    decayHeat(), injectors(), ...), all of which basicNuclearEulerCloud
    inherits unchanged - so the exact same make*Models headers/macros used
    for basicNuclearParcel apply here without modification.

\*---------------------------------------------------------------------------*/

#include "basicNuclearEulerCloud.H"

#include "makeNuclearParcelCloudFunctionObjects.H"

// Kinematic
#include "makeNuclearParcelForces.H" // nuclear variant
#include "makeNuclearEulerParcelForces.H" // two-phase-Euler variant
#include "makeNuclearEulerParcelTurbulenceForces.H" // Brownian motion
#include "makeParcelDispersionModels.H"
#include "makeNuclearEulerParcelChiSwitchModels.H"
#include "makeNuclearParcelInjectionModels.H" // nuclear variant
#include "makeParcelPatchInteractionModels.H"
#include "makeParcelStochasticCollisionModels.H"
#include "makeParcelSurfaceFilmModels.H"

// Thermodynamic-incompressible
#include "makeParcelNuclearHeatTransferModels.H"
#include "makeParcelDecayHeatModels.H"

// MPPIC sub-models
#include "makeMPPICParcelDampingModels.H"
#include "makeMPPICParcelIsotropyModels.H"
#include "makeMPPICParcelPackingModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeNuclearParcelCloudFunctionObjects(basicNuclearEulerCloud);

// Kinematic sub-models
makeNuclearParcelForces(basicNuclearEulerCloud);
makeNuclearEulerParcelForces(basicNuclearEulerCloud);
makeNuclearEulerParcelTurbulenceForces(basicNuclearEulerCloud);
makeParcelDispersionModels(basicNuclearEulerCloud);
makeNuclearEulerParcelChiSwitchModels(basicNuclearEulerCloud);
makeNuclearParcelInjectionModels(basicNuclearEulerCloud);
makeParcelPatchInteractionModels(basicNuclearEulerCloud);
makeParcelStochasticCollisionModels(basicNuclearEulerCloud);
makeParcelSurfaceFilmModels(basicNuclearEulerCloud);

// Nuclear sub-models
makeParcelNuclearHeatTransferModels(basicNuclearEulerCloud);
makeParcelDecayHeatModels(basicNuclearEulerCloud);


// MPPIC sub-models
makeMPPICParcelDampingModels(basicNuclearEulerCloud);
makeMPPICParcelIsotropyModels(basicNuclearEulerCloud);
makeMPPICParcelPackingModels(basicNuclearEulerCloud);


// ************************************************************************* //
