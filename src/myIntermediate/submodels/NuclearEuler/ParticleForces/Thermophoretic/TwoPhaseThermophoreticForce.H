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
    Foam::TwoPhaseThermophoreticForce

Group
    grpLagrangianIntermediateForceSubModels

Description
    Two-phase-Euler-aware particle thermophoretic force.

    Only usable with a cloud whose parcelType::trackingData exposes the
    secondary-phase/chi accessors added by NuclearEulerParcel.

    Applies the same thermophoretic correlation as the plain
    ThermophoreticForce, once against each Euler phase's own temperature
    field (dict entries "T1"/"T2", defaulting to "T.water"/"T.air") and
    carrier thermal conductivity/mean free path (Kc1/lambda1, Kc2/lambda2
    - deliberately separate per phase, since a gas and a liquid carrier
    can have very different molecular mean free paths), and blends the
    two with the particle's own phase-switch parameter chi:
    \f[
        \mathbf{F} = \mathbf{F}_1 (1-\chi) + \mathbf{F}_2 \chi
        \qquad
        \mathbf{F}_i = -D_{t,i}\,\frac{\nabla T_i}{T_{ic}}
    \f]
    The particle's own thermal conductivity Kp is shared between both
    evaluations (it is a property of the particle, not the carrier).
    chi is currently always 0 (no switching model exists yet), so this
    reduces exactly to the primary-phase-only ThermophoreticForce result.

Usage
    \verbatim
    subModels
    {
        particleForces
        {
            twoPhaseThermophoretic;
        }
    }

    twoPhaseThermophoreticCoeffs
    {
        Kp      1.0;
        Kc1     0.6;    // primary (liquid) phase thermal conductivity
        Kc2     0.026;  // secondary (gas) phase thermal conductivity
        lambda1 1e-9;   // primary phase molecular mean free path [m]
        lambda2 6.8e-8; // secondary phase molecular mean free path [m]
        T1      T.water;  // optional, defaults to "T.water"
        T2      T.air;    // optional, defaults to "T.air"
    }
    \endverbatim

SourceFiles
    TwoPhaseThermophoreticForceI.H
    TwoPhaseThermophoreticForce.C

\*---------------------------------------------------------------------------*/

#ifndef TwoPhaseThermophoreticForce_H
#define TwoPhaseThermophoreticForce_H

#include "ParticleForce.H"
#include "volFields.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class TwoPhaseThermophoreticForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class TwoPhaseThermophoreticForce
:
    public ParticleForce<CloudType>
{
protected:

    // Protected data

        //- Name of the primary phase's temperature field
        const word T1Name_;

        //- Name of the secondary phase's temperature field
        const word T2Name_;

        //- Particle thermal conductivity (W/mK)
        const scalar Kp_;

        //- Primary phase thermal conductivity (W/mK)
        const scalar Kc1_;

        //- Secondary phase thermal conductivity (W/mK)
        const scalar Kc2_;

        //- Primary phase molecular mean free path (m)
        const scalar lambda1_;

        //- Secondary phase molecular mean free path (m)
        const scalar lambda2_;

        //- Primary phase temperature gradient interpolator
        autoPtr<interpolation<vector>> gradT1InterpPtr_;

        //- Secondary phase temperature gradient interpolator
        autoPtr<interpolation<vector>> gradT2InterpPtr_;


public:

    //- Runtime type information
    TypeName("twoPhaseThermophoretic");


    // Constructors

        //- Construct from mesh
        TwoPhaseThermophoreticForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& forceType = typeName
        );

        //- Construct copy
        TwoPhaseThermophoreticForce(const TwoPhaseThermophoreticForce& pgf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new TwoPhaseThermophoreticForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~TwoPhaseThermophoreticForce() = default;


    // Member Functions

        // Access

            inline const interpolation<vector>& gradT1Interp() const;

            inline const interpolation<vector>& gradT2Interp() const;


        // Evaluation

            //- Cache fields
            virtual void cacheFields(const bool store);

            //- Calculate the non-coupled force
            virtual forceSuSp calcCoupled
            (
                const typename CloudType::parcelType& p,
                const typename CloudType::parcelType::trackingData& td,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "TwoPhaseThermophoreticForceI.H"

#ifdef NoRepository
    #include "TwoPhaseThermophoreticForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
