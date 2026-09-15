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
    Foam::TwoPhaseGradientDispersionRAS

Description
    Direct port of GradientDispersionRAS to the secondary Euler phase: the
    particle's velocity (as seen against phase 2) is perturbed in the
    direction of -grad(k2), with a Gaussian random number distribution -
    identical formula to the primary-phase model, evaluated with phase 2's
    own k2/epsilon2 fields (see TwoPhaseDispersionRASModel) instead of the
    primary phase's.

\*---------------------------------------------------------------------------*/

#ifndef TwoPhaseGradientDispersionRAS_H
#define TwoPhaseGradientDispersionRAS_H

#include "TwoPhaseDispersionRASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class TwoPhaseGradientDispersionRAS Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class TwoPhaseGradientDispersionRAS
:
    public TwoPhaseDispersionRASModel<CloudType>
{
protected:

    // Protected data

        // Locally cached turbulence fields

            //- Gradient of k2
            const volVectorField* gradk2Ptr_;

            //- Take ownership of grad(k2)
            mutable bool ownGradK2_;


public:

    //- Runtime type information
    TypeName("twoPhaseGradientDispersionRAS");


    // Constructors

        //- Construct from components
        TwoPhaseGradientDispersionRAS
        (
            const dictionary& dict,
            CloudType& owner
        );

        //- Construct copy
        TwoPhaseGradientDispersionRAS
        (
            const TwoPhaseGradientDispersionRAS<CloudType>& dm
        );

        //- Construct and return a clone
        virtual autoPtr<DispersionModel<CloudType>> clone() const
        {
            return autoPtr<DispersionModel<CloudType>>
            (
                new TwoPhaseGradientDispersionRAS<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~TwoPhaseGradientDispersionRAS();


    // Member Functions

        //- Cache carrier fields
        virtual void cacheFields(const bool store);

        //- Update (disperse particles)
        virtual vector update
        (
            const scalar dt,
            const label celli,
            const vector& U,
            const scalar d,
            const scalar rho,
            const vector& Uc,
            const scalar rhoc,
            const scalar muc,
            vector& UTurb,
            scalar& tTurb,
            scalar& tTurbLoc
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TwoPhaseGradientDispersionRAS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
