/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::nuclearParcelInjectionData

Description
    Container class to provide injection data for nuclear parcels

SourceFiles
    nuclearParcelInjectionData.C

\*---------------------------------------------------------------------------*/

#ifndef nuclearParcelInjectionData_H
#define nuclearParcelInjectionData_H

#include "thermoParcelInjectionData.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class nuclearParcelInjectionData;

// Forward declaration of friend functions

Ostream& operator<<
(
    Ostream&,
    const nuclearParcelInjectionData&
);

Istream& operator>>
(
    Istream&,
    nuclearParcelInjectionData&
);

/*---------------------------------------------------------------------------*\
                 Class nuclearParcelInjectionData Declaration
\*---------------------------------------------------------------------------*/

class nuclearParcelInjectionData
:
    public thermoParcelInjectionData
{
//protected:


public:

    //- Runtime type information
    TypeName("nuclearParcelInjectionData");

    // Constructors

        //- Null constructor
        nuclearParcelInjectionData();

        //- Construct from dictionary
        nuclearParcelInjectionData(const dictionary& dict);

        //- Construct from Istream
        nuclearParcelInjectionData(Istream& is);


    //- Destructor
    virtual ~nuclearParcelInjectionData();

    // I-O

        //- Ostream operator
        friend Ostream& operator<<
        (
            Ostream& os,
            const nuclearParcelInjectionData& data
        );

        //- Istream operator
        friend Istream& operator>>
        (
            Istream& is,
            nuclearParcelInjectionData& data
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "nuclearParcelInjectionDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
