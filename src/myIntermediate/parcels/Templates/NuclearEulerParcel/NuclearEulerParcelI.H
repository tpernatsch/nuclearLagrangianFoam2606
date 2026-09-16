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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::NuclearEulerParcel<ParcelType>::constantProperties::
constantProperties()
:
    ParcelType::constantProperties(),
    alphaSmall_(this->dict_, 1e-6)
{}


template<class ParcelType>
inline Foam::NuclearEulerParcel<ParcelType>::constantProperties::
constantProperties
(
    const constantProperties& cp
)
:
    ParcelType::constantProperties(cp),
    alphaSmall_(cp.alphaSmall_)
{}


template<class ParcelType>
inline Foam::NuclearEulerParcel<ParcelType>::constantProperties::
constantProperties
(
    const dictionary& parentDict
)
:
    ParcelType::constantProperties(parentDict),
    alphaSmall_(this->dict_, "alphaSmall", 1e-6)
{}


template<class ParcelType>
inline Foam::NuclearEulerParcel<ParcelType>::NuclearEulerParcel
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    ParcelType(mesh, coordinates, celli, tetFacei, tetPti),
    chi_(0)
{}


template<class ParcelType>
inline Foam::NuclearEulerParcel<ParcelType>::NuclearEulerParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label celli
)
:
    ParcelType(mesh, position, celli),
    chi_(0)
{}


template<class ParcelType>
inline Foam::NuclearEulerParcel<ParcelType>::NuclearEulerParcel
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const scalar dTarget0,
    const vector& U0,
    const vector& f0,
    const vector& angularMomentum0,
    const vector& torque0,
    const constantProperties& constProps
)
:
    ParcelType
    (
        mesh,
        coordinates,
        celli,
        tetFacei,
        tetPti,
        typeId,
        nParticle0,
        d0,
        dTarget0,
        U0,
        f0,
        angularMomentum0,
        torque0,
        constProps
    ),
    chi_(0)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar
Foam::NuclearEulerParcel<ParcelType>::constantProperties::alphaSmall() const
{
    return alphaSmall_.value();
}


// * * * * * * * * * * NuclearEulerParcel Member Functions * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar Foam::NuclearEulerParcel<ParcelType>::chi() const
{
    return chi_;
}


template<class ParcelType>
inline Foam::scalar& Foam::NuclearEulerParcel<ParcelType>::chi()
{
    return chi_;
}


// ************************************************************************* //
