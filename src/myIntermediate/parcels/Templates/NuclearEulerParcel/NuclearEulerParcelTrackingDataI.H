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

template<class ParcelType>
template<class TrackCloudType>
inline Foam::NuclearEulerParcel<ParcelType>::trackingData::trackingData
(
    const TrackCloudType& cloud,
    trackPart part
)
:
    ParcelType::trackingData(cloud, part),
    alpha1Interp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.alpha1()
        )
    ),
    alpha2Interp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.alpha2()
        )
    ),
    U2Interp_
    (
        interpolation<vector>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.U2()
        )
    ),
    rho2Interp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.rho2()
        )
    ),
    mu2Interp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.mu2()
        )
    ),
    T2Interp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.T2()
        )
    ),
    k1Interp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.k1()
        )
    ),
    epsilon1Interp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.epsilon1()
        )
    ),
    k2Interp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.k2()
        )
    ),
    epsilon2Interp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.epsilon2()
        )
    ),
    alpha1c_(Zero),
    alpha2c_(Zero),
    U2c_(Zero),
    rho2c_(Zero),
    mu2c_(Zero),
    T2c_(Zero),
    k1c_(Zero),
    epsilon1c_(Zero),
    k2c_(Zero),
    epsilon2c_(Zero)
{}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::NuclearEulerParcel<ParcelType>::trackingData::alpha1Interp() const
{
    return *alpha1Interp_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::NuclearEulerParcel<ParcelType>::trackingData::alpha2Interp() const
{
    return *alpha2Interp_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::vector>&
Foam::NuclearEulerParcel<ParcelType>::trackingData::U2Interp() const
{
    return *U2Interp_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::NuclearEulerParcel<ParcelType>::trackingData::rho2Interp() const
{
    return *rho2Interp_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::NuclearEulerParcel<ParcelType>::trackingData::mu2Interp() const
{
    return *mu2Interp_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::NuclearEulerParcel<ParcelType>::trackingData::T2Interp() const
{
    return *T2Interp_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::NuclearEulerParcel<ParcelType>::trackingData::k1Interp() const
{
    return *k1Interp_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::NuclearEulerParcel<ParcelType>::trackingData::epsilon1Interp() const
{
    return *epsilon1Interp_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::NuclearEulerParcel<ParcelType>::trackingData::k2Interp() const
{
    return *k2Interp_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::NuclearEulerParcel<ParcelType>::trackingData::epsilon2Interp() const
{
    return *epsilon2Interp_;
}


template<class ParcelType>
inline Foam::scalar
Foam::NuclearEulerParcel<ParcelType>::trackingData::alpha1c() const
{
    return alpha1c_;
}


template<class ParcelType>
inline Foam::scalar&
Foam::NuclearEulerParcel<ParcelType>::trackingData::alpha1c()
{
    return alpha1c_;
}


template<class ParcelType>
inline Foam::scalar
Foam::NuclearEulerParcel<ParcelType>::trackingData::alpha2c() const
{
    return alpha2c_;
}


template<class ParcelType>
inline Foam::scalar&
Foam::NuclearEulerParcel<ParcelType>::trackingData::alpha2c()
{
    return alpha2c_;
}


template<class ParcelType>
inline const Foam::vector&
Foam::NuclearEulerParcel<ParcelType>::trackingData::U2c() const
{
    return U2c_;
}


template<class ParcelType>
inline Foam::vector&
Foam::NuclearEulerParcel<ParcelType>::trackingData::U2c()
{
    return U2c_;
}


template<class ParcelType>
inline Foam::scalar
Foam::NuclearEulerParcel<ParcelType>::trackingData::rho2c() const
{
    return rho2c_;
}


template<class ParcelType>
inline Foam::scalar&
Foam::NuclearEulerParcel<ParcelType>::trackingData::rho2c()
{
    return rho2c_;
}


template<class ParcelType>
inline Foam::scalar
Foam::NuclearEulerParcel<ParcelType>::trackingData::mu2c() const
{
    return mu2c_;
}


template<class ParcelType>
inline Foam::scalar&
Foam::NuclearEulerParcel<ParcelType>::trackingData::mu2c()
{
    return mu2c_;
}


template<class ParcelType>
inline Foam::scalar
Foam::NuclearEulerParcel<ParcelType>::trackingData::T2c() const
{
    return T2c_;
}


template<class ParcelType>
inline Foam::scalar&
Foam::NuclearEulerParcel<ParcelType>::trackingData::T2c()
{
    return T2c_;
}


template<class ParcelType>
inline Foam::scalar
Foam::NuclearEulerParcel<ParcelType>::trackingData::k1c() const
{
    return k1c_;
}


template<class ParcelType>
inline Foam::scalar&
Foam::NuclearEulerParcel<ParcelType>::trackingData::k1c()
{
    return k1c_;
}


template<class ParcelType>
inline Foam::scalar
Foam::NuclearEulerParcel<ParcelType>::trackingData::epsilon1c() const
{
    return epsilon1c_;
}


template<class ParcelType>
inline Foam::scalar&
Foam::NuclearEulerParcel<ParcelType>::trackingData::epsilon1c()
{
    return epsilon1c_;
}


template<class ParcelType>
inline Foam::scalar
Foam::NuclearEulerParcel<ParcelType>::trackingData::k2c() const
{
    return k2c_;
}


template<class ParcelType>
inline Foam::scalar&
Foam::NuclearEulerParcel<ParcelType>::trackingData::k2c()
{
    return k2c_;
}


template<class ParcelType>
inline Foam::scalar
Foam::NuclearEulerParcel<ParcelType>::trackingData::epsilon2c() const
{
    return epsilon2c_;
}


template<class ParcelType>
inline Foam::scalar&
Foam::NuclearEulerParcel<ParcelType>::trackingData::epsilon2c()
{
    return epsilon2c_;
}


// ************************************************************************* //
