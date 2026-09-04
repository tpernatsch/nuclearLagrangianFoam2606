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

#include "albedoFvPatchScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::albedoFvPatchScalarField::albedoFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    diffCoeffName_("diffCoeffAlbedo")
{
    refValue() = 0;
    refGrad() = 0;
    valueFraction() = 1;
}


Foam::albedoFvPatchScalarField::albedoFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    diffCoeffName_(dict.lookupOrDefault<word>("diffCoeffName", "diffCoeffAlbedo"))
{
    beta_ = scalarField("beta", dict, p.size());

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


Foam::albedoFvPatchScalarField::albedoFvPatchScalarField
(
    const albedoFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    beta_(ptf.beta_, mapper),
    diffCoeffName_(ptf.diffCoeffName_)
{
}


Foam::albedoFvPatchScalarField::albedoFvPatchScalarField
(
    const albedoFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    beta_(tppsf.beta_),
    diffCoeffName_(tppsf.diffCoeffName_)
{}


Foam::albedoFvPatchScalarField::albedoFvPatchScalarField
(
    const albedoFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    beta_(tppsf.beta_),
    diffCoeffName_(tppsf.diffCoeffName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::albedoFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);

    beta_.autoMap(m);
}


void Foam::albedoFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const albedoFvPatchScalarField& tiptf =
        refCast<const albedoFvPatchScalarField>(ptf);

    beta_.rmap(tiptf.beta_, addr);
}


void Foam::albedoFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchScalarField& diffCoeffp =
        patch().lookupPatchField<volScalarField, scalar>(diffCoeffName_);

    const scalarField gamma_(0.5*(1 - beta_)/(1 + beta_));

    refGrad() = 0;
    refValue() = 0;

    valueFraction() = gamma_/(gamma_ + diffCoeffp*patch().deltaCoeffs());

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::albedoFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);

    beta_.writeEntry("beta", os);
    os.writeEntryIfDifferent<word>( "diffCoeffName", "diffCoeffAlbedo", diffCoeffName_);

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
        albedoFvPatchScalarField
    );
}

// ************************************************************************* //
