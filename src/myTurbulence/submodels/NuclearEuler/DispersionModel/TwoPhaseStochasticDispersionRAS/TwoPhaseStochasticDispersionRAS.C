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

#include "TwoPhaseStochasticDispersionRAS.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TwoPhaseStochasticDispersionRAS<CloudType>::TwoPhaseStochasticDispersionRAS
(
    const dictionary& dict,
    CloudType& owner
)
:
    TwoPhaseDispersionRASModel<CloudType>(dict, owner)
{}


template<class CloudType>
Foam::TwoPhaseStochasticDispersionRAS<CloudType>::TwoPhaseStochasticDispersionRAS
(
    const TwoPhaseStochasticDispersionRAS<CloudType>& dm
)
:
    TwoPhaseDispersionRASModel<CloudType>(dm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TwoPhaseStochasticDispersionRAS<CloudType>::
~TwoPhaseStochasticDispersionRAS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::vector Foam::TwoPhaseStochasticDispersionRAS<CloudType>::update
(
    const scalar dt,
    const label celli,
    const vector& U,
    const scalar d,
    const scalar rho,
    const vector& Uc,
    const scalar rhoc,
    const scalar muc,
    vector& UTurb,
    scalar& tTurb,
    scalar& tTurbLoc
)
{
    Random& rnd = this->owner().rndGen();

    const scalar cps = 0.16432;

    const scalar k = this->k2Field().primitiveField()[celli];
    const scalar epsilon =
        this->epsilon2Field().primitiveField()[celli] + ROOTVSMALL;

    const scalar UrelMag = mag(U - Uc - UTurb);

    tTurbLoc =
        min(k/epsilon, cps*pow(k, 1.5)/epsilon/(UrelMag + SMALL));


    // Parcel is perturbed by the turbulence
    if (dt < tTurbLoc)
    {
        tTurb += dt;

        if (tTurb > tTurbLoc)
        {
            tTurb = 0;

            const scalar sigma = sqrt(2*k/3.0);

            // Calculate a random direction dir distributed uniformly
            // in spherical coordinates

            const scalar theta = rnd.sample01<scalar>()*twoPi;
            const scalar u = 2*rnd.sample01<scalar>() - 1;

            const scalar a = sqrt(1 - sqr(u));
            const vector dir(a*cos(theta), a*sin(theta), u);

            UTurb = sigma*mag(rnd.GaussNormal<scalar>())*dir;
        }
    }
    else
    {
        tTurb = GREAT;
        UTurb = Zero;
    }

    return Uc + UTurb;
}


// ************************************************************************* //
