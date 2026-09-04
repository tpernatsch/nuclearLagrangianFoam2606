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

#include "depositionFvPatchScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::depositionFvPatchScalarField::depositionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
   diffCoeffName_("D" + this->internalField().name()) 
{
    refValue() = 0;
    refGrad() = 0;
    valueFraction() = 1;
}


Foam::depositionFvPatchScalarField::depositionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
     diffCoeffName_(dict.lookupOrDefault<word>("diffCoeffName", "D" + this->internalField().name()))
{
    gamma_ = scalarField("gamma", dict, p.size());
    c0_ = scalarField("c0", dict, p.size());

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0;
        valueFraction() = 1;
    }
}


Foam::depositionFvPatchScalarField::depositionFvPatchScalarField
(
    const depositionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    gamma_(ptf.gamma_, mapper),
    c0_(ptf.c0_, mapper)
{
}


Foam::depositionFvPatchScalarField::depositionFvPatchScalarField
(
    const depositionFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    gamma_(tppsf.gamma_),
    c0_(tppsf.c0_),
    diffCoeffName_(tppsf.diffCoeffName_)
{}


Foam::depositionFvPatchScalarField::depositionFvPatchScalarField
(
    const depositionFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    gamma_(tppsf.gamma_),
    c0_(tppsf.c0_),
    diffCoeffName_(tppsf.diffCoeffName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::depositionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);

    gamma_.autoMap(m);
    c0_.autoMap(m);
}


void Foam::depositionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const depositionFvPatchScalarField& tiptf =
        refCast<const depositionFvPatchScalarField>(ptf);

    gamma_.rmap(tiptf.gamma_, addr);
    c0_.rmap(tiptf.c0_, addr);
}


void Foam::depositionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchScalarField& diffCoeffp =
        patch().lookupPatchField<volScalarField, scalar>(diffCoeffName_);

    refGrad() = 0.0;
    refValue() = c0_;

    valueFraction() = gamma_/(diffCoeffp*patch().deltaCoeffs());

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::depositionFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);

    gamma_.writeEntry("gamma", os);
    c0_.writeEntry("c0", os);
    os.writeEntryIfDifferent<word>("diffCoeffName", "diffCoeffDeposition", diffCoeffName_);

    refValue().writeEntry("refValue", os);
    refGrad().writeEntry("refGradient", os);
    valueFraction().writeEntry("valueFraction", os);
    fvPatchScalarField::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        depositionFvPatchScalarField
    );
}

// ************************************************************************* //
