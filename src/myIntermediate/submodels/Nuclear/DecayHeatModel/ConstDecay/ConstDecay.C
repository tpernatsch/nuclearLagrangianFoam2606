/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "ConstDecay.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ConstDecay<CloudType>::ConstDecay
(
    const dictionary& dict,
    CloudType& cloud
)
:
    DecayHeatModel<CloudType>(dict, cloud, typeName),
    decayPower_(readScalar(this->coeffDict().lookup("decayPower"))),
    d0_(3.703e-10),
    Q0_(2.146e-21)
{}


template<class CloudType>
Foam::ConstDecay<CloudType>::ConstDecay(const ConstDecay<CloudType>& dhm)
:
    DecayHeatModel<CloudType>(dhm),
    decayPower_(dhm.decayPower_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::ConstDecay<CloudType>::decayPower(const scalar d_part) const
{
    // Calculate decay power: (d_part/d0_)^3 * Q0_/100
    scalar diameterRatio = d_part / d0_;
    return pow(diameterRatio, 3) * Q0_ / 100.0;
}


// ************************************************************************* //
