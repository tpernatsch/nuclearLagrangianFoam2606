/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
inline Foam::NuclearParcel<ParcelType>::trackingData::trackingData
(
    const TrackCloudType& cloud,
    trackPart part
)
:
    ParcelType::trackingData(cloud, part),
    Cp_(cloud.constProps().Cp0()),
    kappa_(cloud.constProps().kappa0()),
    TInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.T()
        )
    ),
    GInterp_(nullptr),
    Tc_(Zero)
{
    if (cloud.radiation())
    {
        GInterp_.reset
        (
            interpolation<scalar>::New
            (
                cloud.solution().interpolationSchemes(),
                cloud.mesh().objectRegistry::template
                    lookupObject<volScalarField>("G")
            ).ptr()
        );
    }
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::NuclearParcel<ParcelType>::trackingData::TInterp() const
{
    return *TInterp_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::NuclearParcel<ParcelType>::trackingData::GInterp() const
{
    if (!GInterp_)
    {
        FatalErrorInFunction
            << "Radiation G interpolation object not set"
            << abort(FatalError);
    }

    return *GInterp_;
}


template<class ParcelType>
inline Foam::scalar Foam::NuclearParcel<ParcelType>::trackingData::Tc() const
{
    return Tc_;
}


template<class ParcelType>
inline Foam::scalar& Foam::NuclearParcel<ParcelType>::trackingData::Tc()
{
    return Tc_;
}


template<class ParcelType>
inline Foam::scalar Foam::NuclearParcel<ParcelType>::trackingData::Cp() const
{
    return Cp_;
}


template<class ParcelType>
inline Foam::scalar Foam::NuclearParcel<ParcelType>::trackingData::kappa() const
{
    return kappa_;
}


// ************************************************************************* //
