/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::CrankNicolson

Description
    Crank-Nicolson integration scheme.

    \f[
        \Delta \phi = (A - B \phi^n) \frac{\Delta t}{1 + \frac{B}{2} \Delta t}
    \f]

\*---------------------------------------------------------------------------*/

#ifndef CrankNicolson_H
#define CrankNicolson_H

#include "integrationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace integrationSchemes
{

/*---------------------------------------------------------------------------*\
                           Class CrankNicolson Declaration
\*---------------------------------------------------------------------------*/

class CrankNicolson
:
    public integrationScheme
{
public:

    //- Runtime type information
    TypeName("CrankNicolson");


    // Constructors

        //- Construct
        CrankNicolson();

        //- Construct and return clone
        virtual autoPtr<integrationScheme> clone() const
        {
            return autoPtr<integrationScheme>(new CrankNicolson(*this));
        }


    //- Destructor
    virtual ~CrankNicolson();


    // Member Functions

        //- Return the integration effective time step
        virtual scalar dtEff(const scalar dt, const scalar Beta) const;

        //- Return the integral of the effective time step (using an CrankNicolson
        //  integration method)
        virtual scalar sumDtEff(const scalar dt, const scalar Beta) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace integrationSchemes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
