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
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::reactorParameters::output
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
)
{
    OFstream& file = this->file();

    file << tab << minValue << tab << maxValue;

    Log << "    " << outputName << " min: " << minValue 
        << " max: " << maxValue << endl;
}


template<class Type>
void Foam::functionObjects::reactorParameters::calcReactorParameters
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const word& outputFieldName
)
{
    const label proci = Pstream::myProcNo();

 List<Type> minVs(Pstream::nProcs(), pTraits<Type>::max);
    labelList minCells(Pstream::nProcs(), Zero);
    List<vector> minCs(Pstream::nProcs(), Zero);

    List<Type> maxVs(Pstream::nProcs(), pTraits<Type>::min);
    labelList maxCells(Pstream::nProcs(), Zero);
    List<vector> maxCs(Pstream::nProcs(), Zero);

    labelPair minMaxIds = findMinMax(field);

    label minId = minMaxIds.first();
    if (minId != -1)
    {
        minVs[proci] = field[minId];
        minCells[proci] = minId;
        minCs[proci] = mesh_.C()[minId];
    }

    label maxId = minMaxIds.second();
    if (maxId != -1)
    {
        maxVs[proci] = field[maxId];
        maxCells[proci] = maxId;
        maxCs[proci] = mesh_.C()[maxId];
    }

// Find min/max boundary field info
    const auto& fieldBoundary = field.boundaryField();
    const auto& CfBoundary = mesh_.C().boundaryField();

    forAll(fieldBoundary, patchi)
    {
        const Field<Type>& fp = fieldBoundary[patchi];
        if (fp.size())
        {
            const vectorField& Cfp = CfBoundary[patchi];

            const labelList& faceCells =
                fieldBoundary[patchi].patch().faceCells();

            minMaxIds = findMinMax(fp);

            minId = minMaxIds.first();
            if (minVs[proci] > fp[minId])
            {
                minVs[proci] = fp[minId];
                minCells[proci] = faceCells[minId];
                minCs[proci] = Cfp[minId];
            }

            maxId = minMaxIds.second();
            if (maxVs[proci] < fp[maxId])
            {
                maxVs[proci] = fp[maxId];
                maxCells[proci] = faceCells[maxId];
                maxCs[proci] = Cfp[maxId];
            }
        }
    }

    // Collect info from all processors and output
    Pstream::allGatherList(minVs);
    Pstream::allGatherList(minCells);
    Pstream::allGatherList(minCs);

    Pstream::allGatherList(maxVs);
    Pstream::allGatherList(maxCells);
    Pstream::allGatherList(maxCs);

    minId = findMin(minVs);
    const Type& minValue = minVs[minId];
    const label minCell = minCells[minId];
    const vector& minC = minCs[minId];

    maxId = findMax(maxVs);
    const Type& maxValue = maxVs[maxId];
    const label maxCell = maxCells[maxId];
    const vector& maxC = maxCs[maxId];

    output
    (
        field.name(),
        outputFieldName,
        minCell,
        maxCell,
        minC,
        maxC,
        minId,
        maxId,
        minValue,
        maxValue
    );
}


template<class Type>
void Foam::functionObjects::reactorParameters::calcReactorParameters
(
    const word& fieldName,
    const modeType& mode
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (obr_.foundObject<fieldType>(fieldName))
    {
        const fieldType& field = lookupObject<fieldType>(fieldName);

        switch (mode)
        {
            case mdMag:
            {
                calcReactorParameters<scalar>
                (
                    mag(field),
                    word("mag(" + fieldName + ")")
                );
                break;
            }
            case mdCmpt:
            {
                calcReactorParameters(field, fieldName);
                break;
            }
            default:
            {
                FatalErrorInFunction
                    << "Unknown mode: " << modeTypeNames_[mode_]
                    << exit(FatalError);
            }
        }
    }
}



// ************************************************************************* //
