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

#include "TwoPhaseSaffmanMeiLiftForce.H"
#include "mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::TwoPhaseSaffmanMeiLiftForce<CloudType>::Cl
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const vector& curlUc,
    const scalar rhoc,
    const scalar Re,
    const scalar muc
) const
{
    scalar Rew = rhoc*mag(curlUc)*sqr(p.d())/(muc + ROOTVSMALL);
    scalar beta = 0.5*(Rew/(Re + ROOTVSMALL));
    scalar alpha = 0.3314*sqrt(beta);
    scalar f = (1.0 - alpha)*exp(-0.1*Re) + alpha;

    scalar Cld = 0.0;
    if (Re < 40)
    {
        Cld = 6.46*f;
    }
    else
    {
        Cld = 6.46*0.0524*sqrt(beta*Re);
    }

    return 3.0/(mathematical::twoPi*sqrt(Rew + ROOTVSMALL))*Cld;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TwoPhaseSaffmanMeiLiftForce<CloudType>::TwoPhaseSaffmanMeiLiftForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& forceType
)
:
    TwoPhaseLiftForce<CloudType>(owner, mesh, dict, forceType)
{}


template<class CloudType>
Foam::TwoPhaseSaffmanMeiLiftForce<CloudType>::TwoPhaseSaffmanMeiLiftForce
(
    const TwoPhaseSaffmanMeiLiftForce& lf
)
:
    TwoPhaseLiftForce<CloudType>(lf)
{}


// ************************************************************************* //
