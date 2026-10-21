/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
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

#include "McLaughlinMeiLiftForce.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::McLaughlinMeiLiftForce<CloudType>::J
(
    const scalar epsilon
)
{
    // Mei (1992), fit of McLaughlin (1991); J_inf = 2.255, 0.3 J_inf = 0.6765.
    // epsilon bounded away from 0 (log10); J -> 0 there anyway
    const scalar eps = max(epsilon, 1e-6);

    return
        0.3*2.255
       *(1.0 + tanh(2.5*(log10(eps) + 0.191)))
       *(2.0/3.0 + tanh(6.0*eps - 1.92));
}


template<class CloudType>
Foam::scalar Foam::McLaughlinMeiLiftForce<CloudType>::Cl
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const vector& curlUc,
    const scalar Re,
    const scalar muc
) const
{
    // shear Reynolds number and McLaughlin's epsilon = sqrt(Re_omega)/Re
    const scalar Rew = td.rhoc()*mag(curlUc)*sqr(p.d())/(muc + ROOTVSMALL);
    const scalar epsilon = sqrt(Rew)/(Re + ROOTVSMALL);

    // Saffman's coefficient 27 J_inf/(2 pi^2 sqrt(Re_omega)) with J_inf -> J
    return
        27.0*J(epsilon)
       /(2.0*sqr(constant::mathematical::pi)*sqrt(Rew + ROOTVSMALL));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::McLaughlinMeiLiftForce<CloudType>::McLaughlinMeiLiftForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& forceType
)
:
    LiftForce<CloudType>(owner, mesh, dict, forceType)
{}


template<class CloudType>
Foam::McLaughlinMeiLiftForce<CloudType>::McLaughlinMeiLiftForce
(
    const McLaughlinMeiLiftForce& lf
)
:
    LiftForce<CloudType>(lf)
{}


// ************************************************************************* //
