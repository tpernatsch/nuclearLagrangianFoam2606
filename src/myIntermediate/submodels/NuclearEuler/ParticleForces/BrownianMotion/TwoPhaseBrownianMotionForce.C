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

#include "TwoPhaseBrownianMotionForce.H"
#include "mathematicalConstants.H"
#include "fundamentalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TwoPhaseBrownianMotionForce<CloudType>::TwoPhaseBrownianMotionForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    lambda_(this->coeffs().getScalar("lambda")),
    turbulence_(this->coeffs().getBool("turbulence"))
{}


template<class CloudType>
Foam::TwoPhaseBrownianMotionForce<CloudType>::TwoPhaseBrownianMotionForce
(
    const TwoPhaseBrownianMotionForce& bmf
)
:
    ParticleForce<CloudType>(bmf),
    lambda_(bmf.lambda_),
    turbulence_(bmf.turbulence_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::TwoPhaseBrownianMotionForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    const scalar dp = p.d();
    const scalar alpha = 2.0*lambda_/dp;
    const scalar cc = 1.0 + alpha*(1.257 + 0.4*exp(-1.1/alpha));

    // Boltzmann constant
    const scalar kb = physicoChemical::k.value();

    Random& rnd = this->owner().rndGen();

    // A single spherical random-direction helper, drawn fresh for each
    // phase's own contribution
    auto sphericalKick = [&](scalar f) -> vector
    {
        const scalar theta = rnd.sample01<scalar>()*constant::mathematical::twoPi;
        const scalar u = 2*rnd.sample01<scalar>() - 1;
        const scalar a = sqrt(max(1 - sqr(u), 0.0));
        const vector dir(a*cos(theta), a*sin(theta), u);

        return f*mag(rnd.GaussNormal<scalar>())*dir;
    };

    // Phase 1 (primary)
    scalar f1 = 0;
    {
        const scalar Tc1 = td.Tc();
        if (turbulence_)
        {
            const scalar Dp1 = kb*Tc1*cc/(3*mathematical::pi*muc*dp);
            f1 = sqrt(2.0*sqr(td.k1c())*sqr(Tc1)/(Dp1*dt + ROOTVSMALL));
        }
        else
        {
            const scalar s0 =
                216*muc*kb*Tc1/(sqr(mathematical::pi)*pow5(dp)*sqr(p.rho())*cc);
            f1 = mass*sqrt(mathematical::pi*s0/dt);
        }
    }

    // Phase 2 (secondary) - same correlation, secondary phase's own
    // temperature/viscosity/turbulence kinetic energy
    scalar f2 = 0;
    {
        const scalar Tc2 = td.T2c();
        const scalar mu2c = td.mu2c();
        if (turbulence_)
        {
            const scalar Dp2 = kb*Tc2*cc/(3*mathematical::pi*mu2c*dp);
            f2 = sqrt(2.0*sqr(td.k2c())*sqr(Tc2)/(Dp2*dt + ROOTVSMALL));
        }
        else
        {
            const scalar s0 =
                216*mu2c*kb*Tc2/(sqr(mathematical::pi)*pow5(dp)*sqr(p.rho())*cc);
            f2 = mass*sqrt(mathematical::pi*s0/dt);
        }
    }

    forceSuSp F1(Zero);
    F1.Su() = sphericalKick(f1);

    forceSuSp F2(Zero);
    F2.Su() = sphericalKick(f2);

    const scalar chi = p.chi();

    return (1.0 - chi)*F1 + chi*F2;
}


// ************************************************************************* //
