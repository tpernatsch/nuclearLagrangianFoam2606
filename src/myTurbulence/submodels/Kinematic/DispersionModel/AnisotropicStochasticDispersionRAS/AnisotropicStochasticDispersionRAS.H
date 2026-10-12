/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::AnisotropicStochasticDispersionRAS

Description
    Discrete-eddy (stochastic) dispersion model, as StochasticDispersionRAS,
    but the eddy velocity fluctuation drawn at each eddy renewal uses an
    anisotropic near-wall RMS distribution below a y+ threshold instead of
    an isotropic sqrt(2k/3) everywhere.

    Below yPlusThreshold, the fluctuation is drawn independently in a local
    wall-oriented frame (streamwise, wall-normal, spanwise) using the
    near-wall RMS correlation of Dreeben and Pope (1997):

        u+_rms = 0.4*y+ / (1 + 0.0239*y+^1.496)
        v+_rms = 0.116*y+^2 / (1 + 0.0203*y+ + 0.0014*y+^2.421)
        w+_rms = 0.19*y+ / (1 + 0.0361*y+^1.322)

    (Dreeben, T.D., Pope, S.B., 1997. Probability density function and
    Reynolds-stress modeling of near-wall turbulent flows. Phys. Fluids 9,
    154-163.)

    The streamwise direction is taken as the local mean flow direction
    projected onto the wall-tangent plane; where that projection is
    degenerate (mean flow ~ wall-normal, e.g. a stagnation point), an
    arbitrary stable tangential direction is used instead.

    Above yPlusThreshold, the draw is unchanged from StochasticDispersionRAS
    (isotropic, sqrt(2k/3)).

    Note: y+ and the friction velocity are computed directly from the
    near-wall velocity gradient (as for a wall-resolved, non-wall-function
    mesh). This model is only meaningful on a mesh that actually resolves
    the near-wall region (y+ ~ O(1) at the first cell) - the correlation
    above is not applicable behind a log-law wall function.

    Does not include the mean-drift correction term used in the
    corresponding continuous (Langevin/CRW) near-wall formulations - left
    for a later revision.

\*---------------------------------------------------------------------------*/

#ifndef AnisotropicStochasticDispersionRAS_H
#define AnisotropicStochasticDispersionRAS_H

#include "DispersionRASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class AnisotropicStochasticDispersionRAS Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class AnisotropicStochasticDispersionRAS
:
    public DispersionRASModel<CloudType>
{
    // Private Data

        //- y+ below which the anisotropic near-wall correlation is used;
        //  above it, the isotropic sqrt(2k/3) draw is used instead
        scalar yPlusThreshold_;

        //- Distance to nearest wall
        autoPtr<volScalarField> yPtr_;

        //- Friction velocity of the nearest wall face, transported into
        //  the volume
        autoPtr<volScalarField> uTauPtr_;

        //- Kinematic viscosity
        autoPtr<volScalarField> nuPtr_;

        //- Unit normal (into the fluid) of the nearest wall face,
        //  transported into the volume
        autoPtr<volVectorField> nWallPtr_;


public:

    //- Runtime type information
    TypeName("anisotropicStochasticDispersionRAS");


    // Constructors

        //- Construct from components
        AnisotropicStochasticDispersionRAS
        (
            const dictionary& dict,
            CloudType& owner
        );

        //- Construct copy
        AnisotropicStochasticDispersionRAS
        (
            const AnisotropicStochasticDispersionRAS<CloudType>& dm
        );

        //- Construct and return a clone
        virtual autoPtr<DispersionModel<CloudType>> clone() const
        {
            return autoPtr<DispersionModel<CloudType>>
            (
                new AnisotropicStochasticDispersionRAS<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~AnisotropicStochasticDispersionRAS();


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
    #include "AnisotropicStochasticDispersionRAS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
