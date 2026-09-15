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
    Foam::TwoPhaseBrownianMotionForce

Group
    grpLagrangianIntermediateForceSubModels

Description
    Two-phase-Euler-aware particle Brownian motion force.

    Only usable with a cloud whose parcelType::trackingData exposes the
    secondary-phase/chi accessors added by NuclearEulerParcel.

    Reference (correlation, applied once per phase):
    \verbatim
        Li, A., & Ahmadi, G. (1992).
        Dispersion and deposition of spherical particles from point sources
        in a turbulent channel flow.
        Aerosol science and technology, 16(4), 209-226.
    \endverbatim

    Unlike myTurbulence's single-phase BrownianMotionForce (which looks
    its turbulence kinetic energy up directly from a registered
    turbulence model), this class reuses the primary/secondary phase
    turbulence kinetic energy already interpolated into trackingData by
    NuclearEulerParcel (td.k1c()/td.k2c()) - no separate mesh lookup is
    needed here.

    Two independent isotropic random kicks are drawn, one per phase
    (each using that phase's own temperature/viscosity/turbulence
    kinetic energy), and blended with the particle's own phase-switch
    parameter chi:
    \f[
        \mathbf{F} = \mathbf{F}_1 (1-\chi) + \mathbf{F}_2 \chi
    \f]
    chi is currently always 0 (no switching model exists yet), so this
    reduces exactly to a single-phase Brownian force evaluated against
    the primary phase only.

    Entries
    \table
        Property   | Description                       | Type   | Reqd | Deflt
        lambda     | Molecular free path length         | scalar | yes  | -
        turbulence | Use turbulence                     | bool   | yes  | -
    \endtable

Usage
    \verbatim
    subModels
    {
        particleForces
        {
            twoPhaseBrownianMotion;
        }
    }

    twoPhaseBrownianMotionCoeffs
    {
        lambda      1e-9;
        turbulence  yes;
    }
    \endverbatim

SourceFiles
    TwoPhaseBrownianMotionForce.C

\*---------------------------------------------------------------------------*/

#ifndef TwoPhaseBrownianMotionForce_H
#define TwoPhaseBrownianMotionForce_H

#include "ParticleForce.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class TwoPhaseBrownianMotionForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class TwoPhaseBrownianMotionForce
:
    public ParticleForce<CloudType>
{
    // Private Data

        //- Molecular free path length [m]
        const scalar lambda_;

        //- Turbulence flag
        const bool turbulence_;


public:

    //- Runtime type information
    TypeName("twoPhaseBrownianMotion");


    // Constructors

        //- Construct from mesh
        TwoPhaseBrownianMotionForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        TwoPhaseBrownianMotionForce(const TwoPhaseBrownianMotionForce& bmf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new TwoPhaseBrownianMotionForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~TwoPhaseBrownianMotionForce() = default;


    // Member Functions

        // Access

            //- The molecular free path length [m]
            scalar lambda() const noexcept { return lambda_; }

            //- The turbulence flag
            bool turbulence() const noexcept { return turbulence_; }


        // Evaluation

            //- Calculate the coupled force
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

#ifdef NoRepository
    #include "TwoPhaseBrownianMotionForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
