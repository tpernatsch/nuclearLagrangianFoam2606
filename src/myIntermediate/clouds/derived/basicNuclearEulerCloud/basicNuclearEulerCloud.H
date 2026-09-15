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
    Foam::basicNuclearEulerCloud

Description
    Cloud class to introduce two-phase-Euler-aware thermodynamic and
    nuclear parcels: NuclearCloud with an additional, secondary Euler
    carrier phase (see NuclearEulerCloud.H).

\*---------------------------------------------------------------------------*/

#ifndef basicNuclearEulerCloud_H
#define basicNuclearEulerCloud_H

#include "NuclearEulerCloud.H"
#include "basicNuclearEulerParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef NuclearEulerCloud
    <
        NuclearCloud
        <
            KinematicCloud
            <
                Cloud
                <
                    basicNuclearEulerParcel
                >
            >
        >
    > basicNuclearEulerCloud;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
