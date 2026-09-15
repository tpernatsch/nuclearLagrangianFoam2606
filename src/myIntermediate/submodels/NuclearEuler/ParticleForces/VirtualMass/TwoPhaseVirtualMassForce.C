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

#include "TwoPhaseVirtualMassForce.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TwoPhaseVirtualMassForce<CloudType>::TwoPhaseVirtualMassForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& forceType
)
:
    TwoPhasePressureGradientForce<CloudType>(owner, mesh, dict, forceType),
    Cvm_(this->coeffs().getScalar("Cvm"))
{}


template<class CloudType>
Foam::TwoPhaseVirtualMassForce<CloudType>::TwoPhaseVirtualMassForce
(
    const TwoPhaseVirtualMassForce& vmf
)
:
    TwoPhasePressureGradientForce<CloudType>(vmf),
    Cvm_(vmf.Cvm_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::TwoPhaseVirtualMassForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value =
        TwoPhasePressureGradientForce<CloudType>::calcCoupled
        (
            p, td, dt, mass, Re, muc
        );

    value.Su() *= Cvm_;

    return value;
}


template<class CloudType>
Foam::scalar Foam::TwoPhaseVirtualMassForce<CloudType>::massAdd
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar mass
) const
{
    const scalar chi = p.chi();

    const scalar massAdd1 = mass*td.rhoc()/p.rho()*Cvm_;
    const scalar massAdd2 = mass*td.rho2c()/p.rho()*Cvm_;

    return (1.0 - chi)*massAdd1 + chi*massAdd2;
}


// ************************************************************************* //
