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

#include "TwoPhaseSphereDragForce.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::TwoPhaseSphereDragForce<CloudType>::CdRe
(
    const scalar Re
) const
{
    // (AOB:Eq. 35)
    if (Re > 1000.0)
    {
        return 0.424*Re;
    }

    return 24.0*(1.0 + (1.0/6.0)*pow(Re, 2.0/3.0));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TwoPhaseSphereDragForce<CloudType>::TwoPhaseSphereDragForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, false)
{}


template<class CloudType>
Foam::TwoPhaseSphereDragForce<CloudType>::TwoPhaseSphereDragForce
(
    const TwoPhaseSphereDragForce<CloudType>& df
)
:
    ParticleForce<CloudType>(df)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::TwoPhaseSphereDragForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    // Phase 1 (primary), identical to plain SphereDragForce
    const scalar Sp1 = mass*0.75*muc*CdRe(Re)/(p.rho()*sqr(p.d()));

    // Phase 2 (secondary), same rho0/mu0 convention, slip vs U2c
    const scalar rho0 = this->owner().constProps().rho0();
    const scalar Re2 = max(p.Re(rho0, p.U(), td.U2c(), p.d(), muc), SMALL);
    const scalar Sp2 = mass*0.75*muc*CdRe(Re2)/(p.rho()*sqr(p.d()));

    const scalar chi = p.chi();

    return forceSuSp(Zero, (1.0 - chi)*Sp1 + chi*Sp2);
}


// ************************************************************************* //
