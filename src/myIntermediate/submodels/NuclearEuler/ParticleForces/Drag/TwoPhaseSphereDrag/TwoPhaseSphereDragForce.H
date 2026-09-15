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
    Foam::TwoPhaseSphereDragForce

Group
    grpLagrangianIntermediateForceSubModels

Description
    Two-phase-Euler-aware sphere-drag force.

    Only usable with a cloud whose parcelType::trackingData exposes the
    secondary-phase/chi accessors added by NuclearEulerParcel.

    Applies the same sphere-drag correlation as the plain SphereDragForce
    (Putnam 1961 / AOB Eq. 34-35), once against each Euler phase, and
    blends the two linearised drag coefficients with the particle's own
    phase-switch parameter chi:
    \f[
        \mathrm{Sp} = \mathrm{Sp}_1 (1-\chi) + \mathrm{Sp}_2 \chi
        \qquad
        \mathrm{Sp}_i =
            \frac{3}{4}\frac{\mu_{ic}\,\mathrm{C}_\mathrm{D}(\mathrm{Re}_i)\,
            \mathrm{Re}_i}{\rho_p \, d_p^2}\, m
    \f]
    Re1/muc are the primary-phase slip Reynolds number/viscosity already
    computed by the caller (identical to plain SphereDragForce); Re2 is
    computed here from the secondary-phase fields cached in trackingData,
    using the same Re() convention NuclearParcel/KinematicParcel already
    use (see the accompanying report).

    chi is currently always 0 (no switching model exists yet), so this
    reduces exactly to the primary-phase-only SphereDragForce result.

Usage
    \verbatim
    subModels
    {
        particleForces
        {
            twoPhaseSphereDrag;
        }
    }
    \endverbatim

SourceFiles
    TwoPhaseSphereDragForce.C

\*---------------------------------------------------------------------------*/

#ifndef TwoPhaseSphereDragForce_H
#define TwoPhaseSphereDragForce_H

#include "ParticleForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class TwoPhaseSphereDragForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class TwoPhaseSphereDragForce
:
    public ParticleForce<CloudType>
{
    // Private Member Functions

        //- Drag coefficient multiplied by Reynolds number
        scalar CdRe(const scalar Re) const;


public:

    //- Runtime type information
    TypeName("twoPhaseSphereDrag");


    // Constructors

        //- Construct from mesh
        TwoPhaseSphereDragForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        TwoPhaseSphereDragForce(const TwoPhaseSphereDragForce<CloudType>& df);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new TwoPhaseSphereDragForce<CloudType>(*this)
            );
        }

        //- No copy assignment
        void operator=(const TwoPhaseSphereDragForce<CloudType>&) = delete;


    //- Destructor
    virtual ~TwoPhaseSphereDragForce() = default;


    // Member Functions

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
    #include "TwoPhaseSphereDragForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
