/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::functionObjects::reactorParameters

Description
    Calculates the value and location of scalar minimum and maximum for a list
    of user-specified fields.

    For variables with a rank greater than zero, either the min/max of a
    component value or the magnitude is reported.  When operating in parallel,
    the processor owning the value is also given.

    Example of function object specification:
    \verbatim
    fieldMinMax1
    {
        type        fieldMinMax;
        libs        ("libfieldFunctionObjects.so");
        ...
        write       yes;
        log         yes;
        location    yes;
        mode        magnitude;
        fields
        (
            U
            p
        );
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: fieldMinMax  | yes         |
        write        | write min/max data to file |  no      | yes
        log          | write min/max data to standard output | no | no
        location     | write location of the min/max value | no | yes
        mode         | calculation mode: magnitude or component | no | magnitude
    \endtable

    Output data is written to the file \<timeDir\>/fieldMinMax.dat

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles

SourceFiles
    reactorParameters.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_reactorParameters_H
#define functionObjects_reactorParameters_H

#include "Switch.H"
#include "Enum.H"
#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "vector.H"
#include "volFieldSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class reactorParameters Declaration
\*---------------------------------------------------------------------------*/

class reactorParameters
:
    public fvMeshFunctionObject,
    public writeFile
{

public:

    enum modeType
    {
        mdMag,
        mdCmpt
    };

protected:

    // Protected data

        //- Mode type names
        static const Enum<modeType> modeTypeNames_;

        //- Mode - only applicable for ranks > 0
        modeType mode_;

        //- Fields to output
         volFieldSelection fieldSet_;


    // Protected Member Functions

        //- Helper function to write the output
        template<class Type>
        void output
        (
            const word& fieldName,
            const word& outputName,
            const label minCell,
            const label maxCell,
            const vector& minC,
            const vector& maxC,
            const label minId,
            const label maxId,
            const Type& minValue,
            const Type& maxValue
        );


	//- Output file header information
        virtual void writeFileHeader(Ostream& os);

        //- Calculate the field
        template<class Type>
        void calcReactorParameters
        (
            const word& fieldName,
            const modeType& mode
        );

        template<class Type>
        void calcReactorParameters
        (
            const GeometricField<Type, fvPatchField, volMesh>& field,
            const word& outputFieldName
        );


public:

    //- Runtime type information
    TypeName("reactorParameters");


    // Constructors

        //- Construct from Time and dictionary
        reactorParameters
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


 //- No copy construct
        reactorParameters(const reactorParameters&) = delete;

        //- No copy assignment
        void operator=(const reactorParameters&) = delete;


    //- Destructor
    virtual ~reactorParameters() = default;


    // Member Functions

        //- Read the field data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the reactorParameters
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "reactorParametersTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
