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

template<class CloudType>
inline const Foam::NuclearEulerCloud<CloudType>&
Foam::NuclearEulerCloud<CloudType>::cloudCopy() const
{
    return cloudCopyPtr_();
}


template<class CloudType>
inline const Foam::volScalarField&
Foam::NuclearEulerCloud<CloudType>::alpha1() const
{
    return alpha1_;
}


template<class CloudType>
inline const Foam::volScalarField&
Foam::NuclearEulerCloud<CloudType>::alpha2() const
{
    return alpha2_;
}


template<class CloudType>
inline const Foam::volVectorField&
Foam::NuclearEulerCloud<CloudType>::U2() const
{
    return U2_;
}


template<class CloudType>
inline const Foam::volScalarField&
Foam::NuclearEulerCloud<CloudType>::rho2() const
{
    return rho2_;
}


template<class CloudType>
inline const Foam::volScalarField&
Foam::NuclearEulerCloud<CloudType>::mu2() const
{
    return mu2_;
}


template<class CloudType>
inline const Foam::volScalarField&
Foam::NuclearEulerCloud<CloudType>::T2() const
{
    return T2_;
}


template<class CloudType>
inline const Foam::volScalarField&
Foam::NuclearEulerCloud<CloudType>::k1() const
{
    return k1_;
}


template<class CloudType>
inline const Foam::volScalarField&
Foam::NuclearEulerCloud<CloudType>::epsilon1() const
{
    return epsilon1_;
}


template<class CloudType>
inline const Foam::volScalarField&
Foam::NuclearEulerCloud<CloudType>::k2() const
{
    return k2_;
}


template<class CloudType>
inline const Foam::volScalarField&
Foam::NuclearEulerCloud<CloudType>::epsilon2() const
{
    return epsilon2_;
}


template<class CloudType>
inline const Foam::volScalarField&
Foam::NuclearEulerCloud<CloudType>::d2() const
{
    return d2_;
}


template<class CloudType>
inline bool Foam::NuclearEulerCloud<CloudType>::k1IsRAS() const
{
    return k1IsRAS_;
}


template<class CloudType>
inline bool Foam::NuclearEulerCloud<CloudType>::k2IsRAS() const
{
    return k2IsRAS_;
}


template<class CloudType>
inline const Foam::DispersionModel<typename CloudType::kinematicCloudType>&
Foam::NuclearEulerCloud<CloudType>::dispersion2() const
{
    return *dispersion2_;
}


template<class CloudType>
inline Foam::DispersionModel<typename CloudType::kinematicCloudType>&
Foam::NuclearEulerCloud<CloudType>::dispersion2()
{
    return *dispersion2_;
}


template<class CloudType>
inline const Foam::ChiSwitchModel<typename Foam::NuclearEulerCloud<CloudType>::nuclearEulerCloudType>&
Foam::NuclearEulerCloud<CloudType>::chiSwitchModel() const
{
    return *chiSwitchModel_;
}


template<class CloudType>
inline Foam::ChiSwitchModel<typename Foam::NuclearEulerCloud<CloudType>::nuclearEulerCloudType>&
Foam::NuclearEulerCloud<CloudType>::chiSwitchModel()
{
    return *chiSwitchModel_;
}


// ************************************************************************* //
