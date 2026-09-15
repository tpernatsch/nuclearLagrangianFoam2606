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
    Registers the two-phase-Euler-aware particle force models -
    additional selectable force types, alongside (not instead of) every
    plain single-phase force already registered via makeNuclearParcelForces
    for a given CloudType. Only usable with a cloud whose parcelType
    exposes the secondary-phase/chi accessors added by NuclearEulerParcel
    (i.e. basicNuclearEulerCloud and its derivatives).

\*---------------------------------------------------------------------------*/

#ifndef makeNuclearEulerParcelForces_H
#define makeNuclearEulerParcelForces_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "TwoPhaseGravityForce.H"
#include "TwoPhaseSphereDragForce.H"
#include "TwoPhaseSaffmanMeiLiftForce.H"
#include "TwoPhasePressureGradientForce.H"
#include "TwoPhaseVirtualMassForce.H"
#include "TwoPhaseThermophoreticForce.H"
#include "TwoPhaseBrownianMotionForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeNuclearEulerParcelForces(CloudType)                              \
                                                                              \
    makeParticleForceModelType(TwoPhaseGravityForce, CloudType);             \
    makeParticleForceModelType(TwoPhaseSphereDragForce, CloudType);          \
    makeParticleForceModelType(TwoPhaseSaffmanMeiLiftForce, CloudType);      \
    makeParticleForceModelType(TwoPhasePressureGradientForce, CloudType);    \
    makeParticleForceModelType(TwoPhaseVirtualMassForce, CloudType);         \
    makeParticleForceModelType(TwoPhaseThermophoreticForce, CloudType);      \
    makeParticleForceModelType(TwoPhaseBrownianMotionForce, CloudType);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
