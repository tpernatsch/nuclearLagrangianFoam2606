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

Class
    Foam::basicNuclearEulerParcel

Description
    Definition of the basic two-phase-Euler-aware nuclear parcel:
    NuclearEulerParcel<NuclearParcel<KinematicParcel<particle>>>

SourceFiles
    defineBasicNuclearEulerParcel.C
    makeBasicNuclearEulerParcelSubmodels.C

\*---------------------------------------------------------------------------*/

#ifndef basicNuclearEulerParcel_H
#define basicNuclearEulerParcel_H

#include "contiguous.H"
#include "particle.H"
#include "KinematicParcel.H"
#include "NuclearParcel.H"
#include "NuclearEulerParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef NuclearEulerParcel<NuclearParcel<KinematicParcel<particle>>>
        basicNuclearEulerParcel;

    //- Contiguous data for basicNuclearEulerParcel
    template<> struct is_contiguous<basicNuclearEulerParcel> : std::true_type {};
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
