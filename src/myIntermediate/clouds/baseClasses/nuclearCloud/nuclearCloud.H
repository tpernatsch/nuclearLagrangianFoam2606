/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::nuclearCloud

Description
    Virtual abstract base class for templated NuclearCloud

SourceFiles
    nuclearCloud.C

\*---------------------------------------------------------------------------*/

#ifndef nuclearCloud_H
#define nuclearCloud_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class nuclearCloud Declaration
\*---------------------------------------------------------------------------*/

class nuclearCloud
{
public:

    //- Runtime type information
    TypeName("nuclearCloud");


    // Constructors

        //- Null constructor
        nuclearCloud() = default;


    //- Destructor
    virtual ~nuclearCloud() = default;


    // Member Functions

        // Access

            // Radiation

                //- Return tmp equivalent particulate emission
                virtual tmp<volScalarField> Ep() const = 0;

                //- Return tmp equivalent particulate absorption
                virtual tmp<volScalarField> ap() const = 0;

                //- Return tmp equivalent particulate scattering factor
                virtual tmp<volScalarField> sigmap() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
