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

#include "TwoPhaseGravityForce.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TwoPhaseGravityForce<CloudType>::TwoPhaseGravityForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, false),
    g_(owner.g().value())
{}


template<class CloudType>
Foam::TwoPhaseGravityForce<CloudType>::TwoPhaseGravityForce
(
    const TwoPhaseGravityForce& gf
)
:
    ParticleForce<CloudType>(gf),
    g_(gf.g_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::TwoPhaseGravityForce<CloudType>::calcNonCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    // F1: buoyancy-corrected gravity against the primary phase
    forceSuSp F1(Zero);
    F1.Su() = mass*g_*(1.0 - td.rhoc()/p.rho());

    // F2: the same, against the secondary phase
    forceSuSp F2(Zero);
    F2.Su() = mass*g_*(1.0 - td.rho2c()/p.rho());

    const scalar chi = p.chi();

    return (1.0 - chi)*F1 + chi*F2;
}


// ************************************************************************* //
