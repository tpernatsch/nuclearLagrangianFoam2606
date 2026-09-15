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
    Foam::nuclearEulerCloud

Description
    Virtual, non-templated marker/base class for the templated
    NuclearEulerCloud.

    Mirrors the role that nuclearCloud plays for NuclearCloud: it exists so
    that non-templated code (fvOptions, function objects, RTTI-based type
    queries) can refer to "a NuclearEulerCloud" without needing to know the
    full CloudType template stack it was instantiated with.

    NuclearEulerCloud<CloudType> extends NuclearCloud<CloudType> (which
    already implements the radiation interface declared by nuclearCloud) with
    awareness of a second, Euler-Euler carrier phase: its own velocity,
    density, dynamic viscosity, temperature, turbulence kinetic energy and
    dissipation rate, plus the local volume fractions (void fractions) of
    both the primary (tracking) phase and this secondary phase. See
    NuclearEulerCloud.H for the full description of how these fields are
    used.

SourceFiles
    nuclearEulerCloud.C

\*---------------------------------------------------------------------------*/

#ifndef nuclearEulerCloud_H
#define nuclearEulerCloud_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class eulerCarrierPhaseFields Declaration
\*---------------------------------------------------------------------------*/

//- Bundle of references to one Euler-Euler carrier phase's fields, as
//- seen by a Lagrangian cloud
class eulerCarrierPhaseFields
{
public:

    //- Volume fraction (void fraction) [-]
    const volScalarField& alpha;

    //- Velocity [m/s]
    const volVectorField& U;

    //- Density [kg/m3]
    const volScalarField& rho;

    //- Dynamic viscosity [Pa.s]
    const volScalarField& mu;

    //- Temperature [K]
    const volScalarField& T;

    //- Turbulence kinetic energy [m2/s2]
    const volScalarField& k;

    //- Turbulence dissipation rate [m2/s3]
    const volScalarField& epsilon;

    //- True if this phase's turbulence model is RAS (gates two-phase
    //- dispersion - see NuclearEulerParcel::calcDispersion())
    const bool isRAS;


    //- Construct from components
    eulerCarrierPhaseFields
    (
        const volScalarField& alpha,
        const volVectorField& U,
        const volScalarField& rho,
        const volScalarField& mu,
        const volScalarField& T,
        const volScalarField& k,
        const volScalarField& epsilon,
        const bool isRAS
    )
    :
        alpha(alpha),
        U(U),
        rho(rho),
        mu(mu),
        T(T),
        k(k),
        epsilon(epsilon),
        isRAS(isRAS)
    {}
};


/*---------------------------------------------------------------------------*\
                     Class nuclearEulerCloud Declaration
\*---------------------------------------------------------------------------*/

class nuclearEulerCloud
{
public:

    //- Runtime type information
    TypeName("nuclearEulerCloud");


    // Constructors

        //- Null constructor
        nuclearEulerCloud() = default;


    //- Destructor
    virtual ~nuclearEulerCloud() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
