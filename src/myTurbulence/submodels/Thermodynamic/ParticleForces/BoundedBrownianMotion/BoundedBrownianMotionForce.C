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

#include "BoundedBrownianMotionForce.H"
#include "mathematicalConstants.H"
#include "fundamentalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BoundedBrownianMotionForce<CloudType>::BoundedBrownianMotionForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    rndGen_(owner.rndGen()),
    lambda_(this->coeffs().getScalar("lambda")),
    useSpherical_(this->coeffs().getOrDefault("spherical", true))
{}


template<class CloudType>
Foam::BoundedBrownianMotionForce<CloudType>::BoundedBrownianMotionForce
(
    const BoundedBrownianMotionForce& bmf
)
:
    ParticleForce<CloudType>(bmf),
    rndGen_(bmf.rndGen_),
    lambda_(bmf.lambda_),
    useSpherical_(bmf.useSpherical_)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BoundedBrownianMotionForce<CloudType>::~BoundedBrownianMotionForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::BoundedBrownianMotionForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(Zero);

    const scalar dp = p.d();
    const scalar Tc = td.Tc();
    const scalar rhop = p.rho();

    const scalar alpha = 2.0*lambda_/dp;
    const scalar cc = 1.0 + alpha*(1.257 + 0.4*exp(-1.1/alpha));

    // Boltzmann constant
    const scalar kb = physicoChemical::k.value();

    // Stokes-Cunningham particle momentum relaxation time
    const scalar tauP = rhop*sqr(dp)*cc/(18.0*muc);

    // Effective time step that the cloud's actual velocity integrator
    // (Euler, CrankNicolson, analytical, ...) applies to a rate with
    // implicit coefficient Beta = 1/tauP. Matching the force amplitude to
    // this - rather than to dt directly, as the stock force does - makes
    // the discrete velocity kick reproduce the correct thermal-fluctuation
    // variance for whichever scheme is selected.
    const scalar dtEff = this->owner().UIntegrator().dtEff(dt, 1.0/tauP);

    // 1 - exp(-2*dt/tauP), evaluated to avoid cancellation error for
    // dt/tauP << 1
    const scalar x = dt/tauP;
    const scalar oneMinusExp2x = (x > SMALL) ? (1.0 - exp(-2.0*x)) : 2.0*x;

    // Amplitude that makes Var[deltaU] match the exact Langevin/
    // Ornstein-Uhlenbeck target (kb*Tc/mass)*oneMinusExp2x for every dt
    const scalar f = sqrt(mass*kb*Tc*oneMinusExp2x)/dtEff;

    Random& rnd = this->owner().rndGen();

    if (useSpherical_)
    {
        // To generate a spherical distribution:
        const scalar theta = rnd.sample01<scalar>()*twoPi;
        const scalar u = 2*rnd.sample01<scalar>() - 1;

        const scalar a = sqrt(1 - sqr(u));
        const vector dir(a*cos(theta), a*sin(theta), u);

        value.Su() = f*mag(rnd.GaussNormal<scalar>())*dir;
    }
    else
    {
        // Generate a cubic distribution (3 independent directions)
        value.Su() = f*rnd.GaussNormal<vector>();
    }

    return value;
}


// ************************************************************************* //
