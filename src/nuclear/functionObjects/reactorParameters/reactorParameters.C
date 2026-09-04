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

\*---------------------------------------------------------------------------*/

#include "reactorParameters.H"
#include "fieldTypes.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(reactorParameters, 0);
    addToRunTimeSelectionTable(functionObject, reactorParameters, dictionary);
}
}

const Foam::Enum
<
    Foam::functionObjects::reactorParameters::modeType
>
Foam::functionObjects::reactorParameters::modeTypeNames_
 ({
     { modeType::mdMag,  "magnitude" },
     { modeType::mdCmpt, "component" },
 });


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::reactorParameters::writeFileHeader(Ostream& os)
{

   if (!fieldSet_.updateSelection())
    {
        return;
    }

    if (writtenHeader_)
    {
        writeBreak(file());
    }
    else
    {
        writeHeader(os, "Reactor Parameters");
    }

    writeCommented(os, "Time");


    for (const word& fieldName : fieldSet_.selectionNames())
        {
            writeTabbed(os, fieldName);
        }

    os << endl;

   writtenHeader_ = true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::reactorParameters::reactorParameters
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    mode_(mdMag),
    fieldSet_(mesh_)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::reactorParameters::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    mode_ = modeTypeNames_.getOrDefault("mode", dict, modeType::mdMag);

    fieldSet_.read(dict);

    return true;
}


bool Foam::functionObjects::reactorParameters::execute()
{
    return true;
}


bool Foam::functionObjects::reactorParameters::write()
{

	writeFileHeader(file());

    if (file().good())
    {
        writeCurrentTime(file());
    }

    Log << type() << " " << name() <<  " write:" << nl;


	 for (const word& fieldName : fieldSet_.selectionNames())
    {
        calcReactorParameters<scalar>(fieldName, mdCmpt);
        calcReactorParameters<vector>(fieldName, mode_);
        calcReactorParameters<sphericalTensor>(fieldName, mode_);
        calcReactorParameters<symmTensor>(fieldName, mode_);
        calcReactorParameters<tensor>(fieldName, mode_);
    }

    if (file().good()) 
    {
        file() << endl;
    }
    Log << endl;

    return true;
}



// ************************************************************************* //
