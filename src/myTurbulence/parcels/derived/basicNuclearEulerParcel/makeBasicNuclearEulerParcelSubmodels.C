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
    Registers, in this order:
      - the ordinary (primary-phase) RAS dispersion models, exactly as
        for basicNuclearParcel - a NuclearEulerCloud still needs these for
        its own inherited, unmodified primary-phase dispersionModel entry;
      - the new secondary-phase (TwoPhase*DispersionRAS) models, only
        selectable through NuclearEulerCloud::dispersion2()'s
        dispersionModel2 entry;
      - the two-phase turbulence-dependent particle force
        (TwoPhaseBrownianMotionForce), alongside the plain two-phase
        forces registered by myIntermediate's makeNuclearEulerParcelForces.

\*---------------------------------------------------------------------------*/

#include "basicNuclearEulerCloud.H"

#include "makeParcelTurbulenceDispersionModels.H"
#include "makeNuclearEulerParcelTurbulenceDispersionModels.H"
#include "makeNuclearEulerParcelTurbulenceForces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeParcelTurbulenceDispersionModels(basicNuclearEulerCloud);
    makeNuclearEulerParcelTurbulenceDispersionModels(basicNuclearEulerCloud);
    makeNuclearEulerParcelTurbulenceForces(basicNuclearEulerCloud);
}


// ************************************************************************* //
