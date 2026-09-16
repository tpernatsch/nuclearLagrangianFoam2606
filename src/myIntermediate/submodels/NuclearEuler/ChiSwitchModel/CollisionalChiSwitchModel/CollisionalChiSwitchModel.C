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

#include "CollisionalChiSwitchModel.H"
#include "mathematicalConstants.H"
#include "physicoChemicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class CloudType>
const Foam::Enum
<
    typename Foam::CollisionalChiSwitchModel<CloudType>::collisionModelType
>
Foam::CollisionalChiSwitchModel<CloudType>::collisionModelTypeNames
{
    { collisionModelType::nonTurbulent, "nonTurbulent" },
    { collisionModelType::turbulent, "turbulent" },
};

template<class CloudType>
const Foam::Enum
<
    typename Foam::CollisionalChiSwitchModel<CloudType>::
        interceptionEfficiencyModelType
>
Foam::CollisionalChiSwitchModel<CloudType>::
    interceptionEfficiencyModelTypeNames
{
    { interceptionEfficiencyModelType::cleanBubble, "cleanBubble" },
    {
        interceptionEfficiencyModelType::fullyContaminatedBubble,
        "fullyContaminatedBubble"
    },
    {
        interceptionEfficiencyModelType::partiallyContaminatedBubble,
        "partiallyContaminatedBubble"
    },
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar
Foam::CollisionalChiSwitchModel<CloudType>::turbulentVelocityFluctuation
(
    const scalar epsilon,
    const scalar nu,
    const scalar d,
    const scalar rho,
    const scalar rhoLiquid
)
{
    // mag() on the density contrast: the bubble is lighter than the
    // liquid, so (rho - rhoLiquid) is negative there, and a fractional
    // power (2/3) of a negative number is not defined - the physically
    // relevant quantity is the magnitude of the density contrast anyway.
    return
        0.4*pow(epsilon, 4.0/9.0)*pow(d, 7.0/9.0)/pow(nu, 1.0/3.0)
       *pow(mag(rho - rhoLiquid)/rhoLiquid, 2.0/3.0);
}


template<class CloudType>
Foam::scalar Foam::CollisionalChiSwitchModel<CloudType>::Z
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td
) const
{
    const scalar dp = p.d();
    const scalar db = td.d2c();
    const scalar Nb = 6.0*td.alpha2c()/(mathematical::pi*pow3(db) + ROOTVSMALL);

    if (collisionModel_ == collisionModelType::nonTurbulent)
    {
        const vector Uslip = td.U2c() - p.U();
        return 0.25*mathematical::pi*sqr(db + dp)*mag(Uslip)*Nb;
    }

    // turbulent
    const scalar rhoLiquid = td.rhoc();
    const scalar nuLiquid = td.muc()/rhoLiquid;
    const scalar epsilonLiquid = td.epsilon1c();
    const scalar rhoP = this->owner().constProps().rho0();
    const scalar rhoB = td.rho2c();

    const scalar Up =
        turbulentVelocityFluctuation(epsilonLiquid, nuLiquid, dp, rhoP, rhoLiquid);
    const scalar Ub =
        turbulentVelocityFluctuation(epsilonLiquid, nuLiquid, db, rhoB, rhoLiquid);

    return 5.0*Nb*sqr(db + dp)*sqrt(sqr(Up) + sqr(Ub));
}


template<class CloudType>
Foam::scalar Foam::CollisionalChiSwitchModel<CloudType>::Ec
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td
) const
{
    const scalar dp = p.d();
    const scalar db = td.d2c();
    const scalar R = dp/db;

    const scalar rhoLiquid = td.rhoc();
    const scalar muLiquid = td.muc();
    const scalar Ub = mag(td.U2c());

    const scalar Reb = rhoLiquid*Ub*db/(muLiquid + ROOTVSMALL);

    scalar EcI = 0;
    switch (interceptionEfficiencyModel_)
    {
        case interceptionEfficiencyModelType::cleanBubble:
        {
            const scalar Reb23 = pow(Reb, 2.0/3.0);
            EcI = R*(15.0 + 3.0*Reb23)/(15.0 + Reb23);
            break;
        }
        case interceptionEfficiencyModelType::fullyContaminatedBubble:
        {
            EcI = 1.5*sqr(R)*(1.0 + pow(Reb, 2.0/3.0)/5.0);
            break;
        }
        case interceptionEfficiencyModelType::partiallyContaminatedBubble:
        {
            EcI = sqr(R)*(1.5 + (4.0/15.0)*pow(Reb, 0.72));
            break;
        }
    }

    scalar EcD = 0;
    if (diffusionalEfficiency_ && Reb >= 10.0)
    {
        const scalar kb = physicoChemical::k.value();
        const scalar Tliquid = td.Tc();

        // Stokes-Einstein diffusion coefficient
        const scalar D = kb*Tliquid/(6.0*mathematical::pi*dp*muLiquid);
        const scalar Pe = Ub*db/(D + ROOTVSMALL);

        EcD = 3.192/sqrt(Pe);
    }

    scalar EcG = 0;
    if (gravitationalEfficiency_)
    {
        const scalar rhoP = this->owner().constProps().rho0();
        const scalar gMag = mag(this->owner().g().value());

        // Stokes settling velocity
        const scalar Vs = (rhoP - rhoLiquid)*gMag*sqr(dp)/(18.0*muLiquid);
        const scalar vs = Vs/(Ub + ROOTVSMALL);

        EcG = vs/(1.0 + vs);
    }

    return EcI + EcD + EcG;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CollisionalChiSwitchModel<CloudType>::CollisionalChiSwitchModel
(
    const dictionary& dict,
    CloudType& owner
)
:
    ChiSwitchModel<CloudType>(dict, owner, typeName),
    collisionModel_
    (
        collisionModelTypeNames.get("collisionModel", this->coeffDict())
    ),
    interceptionEfficiencyModel_
    (
        interceptionEfficiencyModelTypeNames.get
        (
            "interceptionEfficiencyModel",
            this->coeffDict()
        )
    ),
    diffusionalEfficiency_(this->coeffDict().getBool("diffusionalEfficiency")),
    gravitationalEfficiency_
    (
        this->coeffDict().getBool("gravitationalEfficiency")
    ),
    Ea_(this->coeffDict().getScalar("Ea"))
{}


template<class CloudType>
Foam::CollisionalChiSwitchModel<CloudType>::CollisionalChiSwitchModel
(
    const CollisionalChiSwitchModel<CloudType>& csm
)
:
    ChiSwitchModel<CloudType>(csm),
    collisionModel_(csm.collisionModel_),
    interceptionEfficiencyModel_(csm.interceptionEfficiencyModel_),
    diffusionalEfficiency_(csm.diffusionalEfficiency_),
    gravitationalEfficiency_(csm.gravitationalEfficiency_),
    Ea_(csm.Ea_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::CollisionalChiSwitchModel<CloudType>::gamma
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td
) const
{
    return Z(p, td)*Ec(p, td)*Ea_;
}


// ************************************************************************* //
