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
    Foam::TwoPhaseGravityForce

Group
    grpLagrangianIntermediateForceSubModels

Description
    Two-phase-Euler-aware particle gravity (buoyancy-corrected) force.

    Only usable with a cloud whose parcelType::trackingData exposes the
    secondary-phase/chi accessors added by NuclearEulerParcel (i.e.
    basicNuclearEulerCloud and its derivatives) - the template only
    compiles once actually instantiated against such a parcel type.

    Computes the same buoyancy-corrected gravity force as the plain
    GravityForce, once against each Euler phase, and blends the two with
    the particle's own phase-switch parameter chi:
    \f[
        \mathbf{F} = \mathbf{F}_1 (1-\chi) + \mathbf{F}_2 \chi
    \f]
    \f[
        \mathbf{F}_1 = m\,\mathbf{g}\left(1 - \frac{\rho_{1c}}{\rho_p}\right)
        \qquad
        \mathbf{F}_2 = m\,\mathbf{g}\left(1 - \frac{\rho_{2c}}{\rho_p}\right)
    \f]
    chi is currently always 0 (no switching model exists yet - see
    NuclearEulerParcel), so this reduces exactly to the primary-phase-only
    GravityForce result until a future switching model sets chi non-zero.

Usage
    \verbatim
    subModels
    {
        particleForces
        {
            twoPhaseGravity;
        }
    }
    \endverbatim

SourceFiles
    TwoPhaseGravityForce.C

\*---------------------------------------------------------------------------*/

#ifndef TwoPhaseGravityForce_H
#define TwoPhaseGravityForce_H

#include "ParticleForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class TwoPhaseGravityForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class TwoPhaseGravityForce
:
    public ParticleForce<CloudType>
{
    // Private Data

        //- Reference to the acceleration due to gravity
        const vector& g_;


public:

    //- Runtime type information
    TypeName("twoPhaseGravity");


    // Constructors

        //- Construct from mesh
        TwoPhaseGravityForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        TwoPhaseGravityForce(const TwoPhaseGravityForce& gf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new TwoPhaseGravityForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~TwoPhaseGravityForce() = default;


    // Member Functions

        // Access

            //- Return the acceleration due to gravity
            const vector& g() const noexcept { return g_; }


        // Evaluation

            //- Calculate the non-coupled force
            virtual forceSuSp calcNonCoupled
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
    #include "TwoPhaseGravityForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
