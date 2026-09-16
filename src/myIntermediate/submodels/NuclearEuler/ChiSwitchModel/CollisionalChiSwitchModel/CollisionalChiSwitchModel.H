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

Class
    Foam::CollisionalChiSwitchModel

Description
    Bubble-particle collision/attachment chi-switch model:
    \f[
        \gamma = Z \, E_c \, E_a
    \f]
    where Z is a bubble-particle collision frequency, E_c a collision
    efficiency (itself a sum of up to three contributions) and E_a a
    (currently constant-only) attachment efficiency.

    Z: one of two mutually exclusive correlations, selected by
    collisionModel:
    \f[
        Z_{nonTurbulent} = \frac{\pi}{4}(d_b+d_p)^2\,|\mathbf{U}_b-\mathbf{U}_p|\,N_b
        \qquad
        Z_{turbulent} = 5 N_b (d_b+d_p)^2 \sqrt{U_p'^2+U_b'^2}
    \f]
    with the turbulent velocity fluctuation of a species i (particle or
    bubble) of diameter d_i and density rho_i, in a liquid of turbulence
    dissipation epsilon and kinematic viscosity nu:
    \f[
        U_i' = \frac{0.4\,\epsilon^{4/9} d_i^{7/9}}{\nu^{1/3}}
               \left(\frac{\rho_i-\rho_{liq}}{\rho_{liq}}\right)^{2/3}
    \f]
    and the bubble number density Nb = 6*alpha2/(pi*d_b^3) (monodisperse
    void-fraction-consistent estimate - a placeholder pending a better
    source for Nb).

    E_c = E_c,i + E_c,D + E_c,g:
      - E_c,i: exactly one of three mutually exclusive correlations,
        selected by interceptionEfficiencyModel (cleanBubble /
        fullyContaminatedBubble / partiallyContaminatedBubble).
      - E_c,D: diffusional contribution, only added when
        diffusionalEfficiency is on AND Re_b >= 10:
        E_c,D = 3.192*Pe^-1/2, Pe = U_b*d_b/D (Stokes-Einstein D).
      - E_c,g: gravitational settling contribution, added when
        gravitationalEfficiency is on: E_c,g = v_s/(1+v_s),
        v_s = V_s/U_b (Stokes settling velocity V_s).

    E_a: currently a single, constant, dictionary-supplied value - a
    placeholder for a future pluggable attachment-efficiency model.

    Re_b = rho_liquid*|U_b|*d_b/mu_liquid throughout.

Usage
    \verbatim
    subModels
    {
        chiSwitchModel  collisional;

        collisionalCoeffs
        {
            collisionModel               turbulent;  // or nonTurbulent
            interceptionEfficiencyModel   cleanBubble; // or
                                          // fullyContaminatedBubble /
                                          // partiallyContaminatedBubble
            diffusionalEfficiency         yes;
            gravitationalEfficiency       yes;
            Ea                            0.5;
        }
    }
    \endverbatim

SourceFiles
    CollisionalChiSwitchModel.C

\*---------------------------------------------------------------------------*/

#ifndef CollisionalChiSwitchModel_H
#define CollisionalChiSwitchModel_H

#include "ChiSwitchModel.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class CollisionalChiSwitchModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CollisionalChiSwitchModel
:
    public ChiSwitchModel<CloudType>
{
public:

    //- Collision frequency (Z) correlation
    enum class collisionModelType
    {
        nonTurbulent,
        turbulent
    };

    //- Collision efficiency (Ec,i) correlation
    enum class interceptionEfficiencyModelType
    {
        cleanBubble,
        fullyContaminatedBubble,
        partiallyContaminatedBubble
    };

    static const Enum<collisionModelType> collisionModelTypeNames;
    static const Enum<interceptionEfficiencyModelType>
        interceptionEfficiencyModelTypeNames;


private:

    // Private Data

        //- Collision frequency correlation
        const collisionModelType collisionModel_;

        //- Collision efficiency (interception) correlation
        const interceptionEfficiencyModelType interceptionEfficiencyModel_;

        //- Include the diffusional collision-efficiency contribution
        //- (only when Re_b >= 10)
        const bool diffusionalEfficiency_;

        //- Include the gravitational collision-efficiency contribution
        const bool gravitationalEfficiency_;

        //- Constant attachment efficiency
        const scalar Ea_;


    // Private Member Functions

        //- Turbulent velocity fluctuation of a species of diameter d and
        //- density rho, suspended in a liquid of dissipation rate
        //- epsilon and kinematic viscosity nu
        static scalar turbulentVelocityFluctuation
        (
            const scalar epsilon,
            const scalar nu,
            const scalar d,
            const scalar rho,
            const scalar rhoLiquid
        );

        //- Bubble-particle collision frequency Z [1/s]
        scalar Z
        (
            const typename CloudType::parcelType& p,
            const typename CloudType::parcelType::trackingData& td
        ) const;

        //- Collision efficiency Ec = Ec,i + Ec,D + Ec,g
        scalar Ec
        (
            const typename CloudType::parcelType& p,
            const typename CloudType::parcelType::trackingData& td
        ) const;


public:

    //- Runtime type information
    TypeName("collisional");


    // Constructors

        //- Construct from components
        CollisionalChiSwitchModel(const dictionary& dict, CloudType& owner);

        //- Construct copy
        CollisionalChiSwitchModel
        (
            const CollisionalChiSwitchModel<CloudType>& csm
        );

        //- Construct and return a clone
        virtual autoPtr<ChiSwitchModel<CloudType>> clone() const
        {
            return autoPtr<ChiSwitchModel<CloudType>>
            (
                new CollisionalChiSwitchModel<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~CollisionalChiSwitchModel() = default;


    // Member Functions

        //- Return the switch rate constant gamma = Z*Ec*Ea [1/s]
        virtual scalar gamma
        (
            const typename CloudType::parcelType& p,
            const typename CloudType::parcelType::trackingData& td
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CollisionalChiSwitchModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
