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
    Foam::TwoPhaseStochasticDispersionRAS

Description
    Direct port of StochasticDispersionRAS to the secondary Euler phase:
    an isotropic random perturbation with variance set by phase 2's own
    k2/epsilon2 fields (see TwoPhaseDispersionRASModel) instead of the
    primary phase's.

\*---------------------------------------------------------------------------*/

#ifndef TwoPhaseStochasticDispersionRAS_H
#define TwoPhaseStochasticDispersionRAS_H

#include "TwoPhaseDispersionRASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class TwoPhaseStochasticDispersionRAS Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class TwoPhaseStochasticDispersionRAS
:
    public TwoPhaseDispersionRASModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("twoPhaseStochasticDispersionRAS");


    // Constructors

        //- Construct from components
        TwoPhaseStochasticDispersionRAS
        (
            const dictionary& dict,
            CloudType& owner
        );

        //- Construct copy
        TwoPhaseStochasticDispersionRAS
        (
            const TwoPhaseStochasticDispersionRAS<CloudType>& dm
        );

        //- Construct and return a clone
        virtual autoPtr<DispersionModel<CloudType>> clone() const
        {
            return autoPtr<DispersionModel<CloudType>>
            (
                new TwoPhaseStochasticDispersionRAS<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~TwoPhaseStochasticDispersionRAS();


    // Member Functions

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
    #include "TwoPhaseStochasticDispersionRAS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
