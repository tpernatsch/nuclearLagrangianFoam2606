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
    Foam::TwoPhasePressureGradientForce

Group
    grpLagrangianIntermediateForceSubModels

Description
    Two-phase-Euler-aware particle pressure-gradient force.

    Only usable with a cloud whose parcelType::trackingData exposes the
    secondary-phase/chi accessors added by NuclearEulerParcel.

    The plain PressureGradientForce computes D(Uc)/Dt = ddt(Uc)+Uc&grad(Uc)
    for one named velocity field "U" (dict entry, default "U") and caches
    it as a registered mesh field/interpolator. This class does the exact
    same thing twice - once for the primary phase's real velocity field
    (dict entry "U1", e.g. "U.water") and once for the secondary phase's
    ("U2", e.g. "U.air") - and blends the two results with the particle's
    own phase-switch parameter chi:
    \f[
        \mathbf{F} = \mathbf{F}_1 (1-\chi) + \mathbf{F}_2 \chi
        \qquad
        \mathbf{F}_i = m\,\frac{\rho_{ic}}{\rho_p}\,\frac{D\mathbf{U}_{ic}}{Dt}
    \f]
    chi is currently always 0 (no switching model exists yet), so this
    reduces exactly to the primary-phase-only PressureGradientForce
    result.

Usage
    \verbatim
    subModels
    {
        particleForces
        {
            twoPhasePressureGradient;
        }
    }

    twoPhasePressureGradientCoeffs
    {
        U1  U.water;   // optional, defaults to "U.water"
        U2  U.air;     // optional, defaults to "U.air"
    }
    \endverbatim

SourceFiles
    TwoPhasePressureGradientForceI.H
    TwoPhasePressureGradientForce.C

\*---------------------------------------------------------------------------*/

#ifndef TwoPhasePressureGradientForce_H
#define TwoPhasePressureGradientForce_H

#include "ParticleForce.H"
#include "volFields.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class TwoPhasePressureGradientForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class TwoPhasePressureGradientForce
:
    public ParticleForce<CloudType>
{
protected:

    // Protected data

        //- Name of the primary phase's velocity field
        const word U1Name_;

        //- Name of the secondary phase's velocity field
        const word U2Name_;

        //- Rate of change of primary phase velocity interpolator
        autoPtr<interpolation<vector>> DU1cDtInterpPtr_;

        //- Rate of change of secondary phase velocity interpolator
        autoPtr<interpolation<vector>> DU2cDtInterpPtr_;


public:

    //- Runtime type information
    TypeName("twoPhasePressureGradient");


    // Constructors

        //- Construct from mesh
        TwoPhasePressureGradientForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& forceType = typeName
        );

        //- Construct copy
        TwoPhasePressureGradientForce(const TwoPhasePressureGradientForce& pgf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new TwoPhasePressureGradientForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~TwoPhasePressureGradientForce() = default;


    // Member Functions

        // Access

            //- Return the primary phase D/Dt(U) interpolator
            inline const interpolation<vector>& DU1cDtInterp() const;

            //- Return the secondary phase D/Dt(U) interpolator
            inline const interpolation<vector>& DU2cDtInterp() const;


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

            //- Return the added mass
            virtual scalar massAdd
            (
                const typename CloudType::parcelType& p,
                const typename CloudType::parcelType::trackingData& td,
                const scalar mass
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "TwoPhasePressureGradientForceI.H"

#ifdef NoRepository
    #include "TwoPhasePressureGradientForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
