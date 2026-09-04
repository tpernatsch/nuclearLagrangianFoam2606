/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "basicNuclearCloud.H"

#include "makeNuclearParcelCloudFunctionObjects.H"

// Kinematic
#include "makeNuclearParcelForces.H" // nuclear variant
#include "makeParcelDispersionModels.H"
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

makeNuclearParcelCloudFunctionObjects(basicNuclearCloud);

// Kinematic sub-models
makeNuclearParcelForces(basicNuclearCloud);
makeParcelDispersionModels(basicNuclearCloud);
makeNuclearParcelInjectionModels(basicNuclearCloud);
makeParcelPatchInteractionModels(basicNuclearCloud);
makeParcelStochasticCollisionModels(basicNuclearCloud);
makeParcelSurfaceFilmModels(basicNuclearCloud);

// Nuclear sub-models
makeParcelNuclearHeatTransferModels(basicNuclearCloud);
makeParcelDecayHeatModels(basicNuclearCloud);


// MPPIC sub-models
makeMPPICParcelDampingModels(basicNuclearCloud);
makeMPPICParcelIsotropyModels(basicNuclearCloud);
makeMPPICParcelPackingModels(basicNuclearCloud);


// ************************************************************************* //
