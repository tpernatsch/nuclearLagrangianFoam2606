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
    Foam::BoundedBrownianMotionForce

Description
    Particle Brownian motion force with a force amplitude that reproduces
    the exact Ornstein-Uhlenbeck velocity-transition variance for whichever
    integrationScheme (Euler, CrankNicolson, analytical, ...) is selected
    for the cloud's velocity update, for every time step size.

    This is a separate force from the stock BrownianMotionForce (Li &
    Ahmadi, 1992), not a replacement for it: the stock amplitude scales as
    dt^(-1/2) unconditionally, which is only valid while dt << tauP (the
    particle momentum relaxation time). For nanoparticles, where
    dt >> tauP is essentially always true, that amplitude diverges as any
    adaptive sub-stepping shrinks dt, feeding back destructively into the
    sub-stepping itself. This force instead matches the discrete velocity
    kick's variance to the exact Langevin/Ornstein-Uhlenbeck target
    (kb*Tc/mp)*(1 - exp(-2*dt/tauP)) using the cloud's own UIntegrator()
    effective time step, so the amplitude stays bounded (and becomes
    dt-independent) in the dt >> tauP regime instead of diverging, while
    reducing to the stock formula in the dt << tauP regime.

    Only the non-turbulence (molecular) formulation is implemented.

    Entries
    \table
        Property | Description                            | Type | Reqd | Deflt
        lambda   | Molecular free path length             | scalar | yes  | -
        spherical | Spherical instead of cubic distribution | bool | no  | true
    \endtable

SourceFiles
    BoundedBrownianMotionForce.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_BoundedBrownianMotionForce_H
#define Foam_BoundedBrownianMotionForce_H

#include "ParticleForce.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class BoundedBrownianMotionForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class BoundedBrownianMotionForce
:
    public ParticleForce<CloudType>
{
    // Private Data

        //- Reference to the cloud random number generator
        Random& rndGen_;

        //- Molecular free path length [m]
        const scalar lambda_;

        //- Spherical vs cubic (transitional)
        bool useSpherical_;


public:

    //- Runtime type information
    TypeName("BoundedBrownianMotion");


    // Constructors

        //- Construct from mesh
        BoundedBrownianMotionForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        BoundedBrownianMotionForce(const BoundedBrownianMotionForce& bmf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new BoundedBrownianMotionForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~BoundedBrownianMotionForce();


    // Member Functions

    // Access

        //- The molecular free path length [m]
        scalar lambda() const noexcept { return lambda_; }


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
    #include "BoundedBrownianMotionForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
