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
    Foam::TwoPhaseLiftForce

Group
    grpLagrangianIntermediateForceSubModels

Description
    Base class for two-phase-Euler-aware particle lift force models.

    Only usable with a cloud whose parcelType::trackingData exposes the
    secondary-phase/chi accessors added by NuclearEulerParcel.

    Mirrors the plain LiftForce, but caches curl(U) for *both* Euler
    phases' real velocity fields (dict entries "U1"/"U2", defaulting to
    "U.water"/"U.air"), computes the lift coefficient Cl() once per phase
    (via the same virtual Cl() hook plain LiftForce uses - concrete
    correlations, e.g. TwoPhaseSaffmanMeiLiftForce, override it exactly
    as before), and blends the two resulting forces with the particle's
    own phase-switch parameter chi:
    \f[
        \mathbf{F} = \mathbf{F}_1 (1-\chi) + \mathbf{F}_2 \chi
        \qquad
        \mathbf{F}_i =
            \frac{m}{\rho_p}\rho_{ic}\,C_{l,i}\,
            \left[(\mathbf{U}_{ic}-\mathbf{U}_p)\times
            (\nabla\times\mathbf{U}_{ic})\right]
    \f]
    chi is currently always 0 (no switching model exists yet), so this
    reduces exactly to the primary-phase-only LiftForce result.

SourceFiles
    TwoPhaseLiftForceI.H
    TwoPhaseLiftForce.C

\*---------------------------------------------------------------------------*/

#ifndef TwoPhaseLiftForce_H
#define TwoPhaseLiftForce_H

#include "ParticleForce.H"
#include "volFields.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class TwoPhaseLiftForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class TwoPhaseLiftForce
:
    public ParticleForce<CloudType>
{
protected:

    // Protected data

        //- Name of the primary phase's velocity field
        const word U1Name_;

        //- Name of the secondary phase's velocity field
        const word U2Name_;

        //- Curl of the primary phase velocity interpolator
        autoPtr<interpolation<vector>> curlU1cInterpPtr_;

        //- Curl of the secondary phase velocity interpolator
        autoPtr<interpolation<vector>> curlU2cInterpPtr_;


    // Protected Member Functions

        //- Calculate the lift coefficient
        virtual scalar Cl
        (
            const typename CloudType::parcelType& p,
            const typename CloudType::parcelType::trackingData& td,
            const vector& curlUc,
            const scalar rhoc,
            const scalar Re,
            const scalar muc
        ) const;


public:

    // Constructors

        //- Construct from mesh
        TwoPhaseLiftForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& forceType
        );

        //- Construct copy
        TwoPhaseLiftForce(const TwoPhaseLiftForce& lf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new TwoPhaseLiftForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~TwoPhaseLiftForce() = default;


    // Member Functions

        // Access

            //- Return the curl of the primary phase velocity interpolator
            inline const interpolation<vector>& curlU1cInterp() const;

            //- Return the curl of the secondary phase velocity
            //- interpolator
            inline const interpolation<vector>& curlU2cInterp() const;


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

#include "TwoPhaseLiftForceI.H"

#ifdef NoRepository
    #include "TwoPhaseLiftForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
