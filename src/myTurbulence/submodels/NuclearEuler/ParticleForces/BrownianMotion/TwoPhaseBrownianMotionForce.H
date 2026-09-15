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

Description
    Two-phase-Euler-aware particle Brownian motion force. Direct port of
    BrownianMotionForce (myTurbulence) to the secondary Euler phase: the
    same correlation is evaluated once per phase, each phase's own
    turbulence kinetic energy located via its own registered turbulence
    model (kModel()/k2Model()), exactly the way the single-phase force
    locates its one turbulence model, and the two results are blended by
    the particle's own phase-switch parameter chi.

    Reference:
    \verbatim
        Li, A., & Ahmadi, G. (1992).
        Dispersion and deposition of spherical particles from point sources
        in a turbulent channel flow.
        Aerosol science and technology, 16(4), 209-226.
    \endverbatim

    Entries
    \table
        Property   | Description                       | Type   | Reqd | Deflt
        lambda     | Molecular free path length         | scalar | yes  | -
        turbulence | Use turbulence                     | bool   | yes  | -
        spherical  | Spherical instead of cubic distribution | bool | no | true
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

        //- Reference to the cloud random number generator
        Random& rndGen_;

        //- Molecular free path length [m]
        const scalar lambda_;

        //- Pointer to the primary phase turbulence kinetic energy field
        const volScalarField* kPtr_;

        //- Pointer to the secondary phase turbulence kinetic energy field
        const volScalarField* k2Ptr_;

        //- Turbulence flag
        const bool turbulence_;

        //- Local ownership of the primary phase turbulence k field
        bool ownK_;

        //- Local ownership of the secondary phase turbulence k field
        bool ownK2_;

        //- Spherical vs cubic (transitional)
        const bool useSpherical_;


    // Private Member Functions

        //- Return the primary phase k field from its turbulence model
        tmp<volScalarField> kModel() const;

        //- Return the secondary phase k field from its turbulence model
        tmp<volScalarField> k2Model() const;


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
    virtual ~TwoPhaseBrownianMotionForce();


    // Member Functions

    // Access

        //- The molecular free path length [m]
        scalar lambda() const noexcept { return lambda_; }

        //- The turbulence flag
        bool turbulence() const noexcept { return turbulence_; }


        // Evaluation

            //- Cache fields
            virtual void cacheFields(const bool store);

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
