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

#include "albedoSP3FvPatchField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::albedoSP3FvPatchField<Type>::albedoSP3FvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    gamma_(scalar(0.5)),
    diffCoeffName_("diffCoeffName"),
    fluxStarAlbedo_("fluxStarAlbedo"),
    forSecondMoment_(false)
{}


    template<class Type>
    Foam::albedoSP3FvPatchField<Type>::albedoSP3FvPatchField
    (
        const fvPatch& p,
        const DimensionedField<Type, volMesh>& iF,
        const dictionary& dict
    )
    :
        fvPatchField<Type>(p, iF, dict),
        gamma_(dict.lookupOrDefault<scalar>("gamma", scalar(0.5))),
        diffCoeffName_(dict.lookupOrDefault<word>("diffCoeffName", "diffCoeffName")),
        fluxStarAlbedo_(dict.lookupOrDefault<word>("fluxStarAlbedo", "fluxStarAlbedo")),
        forSecondMoment_(dict.lookupOrDefault<bool>("forSecondMoment", false))
    {
        evaluate();
    }



template<class Type>
Foam::albedoSP3FvPatchField<Type>::albedoSP3FvPatchField
(
    const albedoSP3FvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
    gamma_(ptf.gamma_),
    diffCoeffName_(ptf.diffCoeffName_),
    fluxStarAlbedo_(ptf.fluxStarAlbedo_),
    forSecondMoment_(ptf.forSecondMoment_)
{
    if (notNull(iF) && mapper.hasUnmapped())
    {
        WarningInFunction
            << "On field " << iF.name() << " patch " << p.name()
            << " patchField " << this->type()
            << " : mapper does not map all values." << nl
            << "    To avoid this warning fully specify the mapping in derived"
            << " patch fields." << endl;
    }
}


template<class Type>
Foam::albedoSP3FvPatchField<Type>::albedoSP3FvPatchField
(
    const albedoSP3FvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),
    gamma_(ptf.gamma_),
    diffCoeffName_(ptf.diffCoeffName_),
    fluxStarAlbedo_(ptf.fluxStarAlbedo_),
    forSecondMoment_(ptf.forSecondMoment_)
{}


template<class Type>
Foam::albedoSP3FvPatchField<Type>::albedoSP3FvPatchField
(
    const albedoSP3FvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    gamma_(ptf.gamma_),
    diffCoeffName_(ptf.diffCoeffName_),
    fluxStarAlbedo_(ptf.fluxStarAlbedo_),
    forSecondMoment_(ptf.forSecondMoment_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::albedoSP3FvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<Type>::autoMap(m);
}


template<class Type>
void Foam::albedoSP3FvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    const albedoSP3FvPatchField<Type>& fgptf =
        refCast<const albedoSP3FvPatchField<Type>>(ptf);
}


template<class Type>
void Foam::albedoSP3FvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const Foam::Field<scalar>& diffCoeff =
        this->patch().template lookupPatchField<volScalarField, scalar>
        (
            diffCoeffName_
        );
    const Foam::Field<scalar> fluxStarAlbedo =
        this->patch().template lookupPatchField<volScalarField, scalar>
        (
            fluxStarAlbedo_
        );

    if(forSecondMoment_)
    {
	    Foam::Field<Type>::operator=
	    (
	        this->patchInternalField()
		+ (this->patchInternalField()*gamma_/diffCoeff*21.0/20.0)/this->patch().deltaCoeffs()
		-  (Type(pTraits<Type>::one)*fluxStarAlbedo*gamma_/diffCoeff*3.0/20.0)/this->patch().deltaCoeffs()
	    );
    }
    else
    {
	    Foam::Field<Type>::operator=
	    (
	        this->patchInternalField()
                + (this->patchInternalField()*gamma_/diffCoeff)/this->patch().deltaCoeffs()
                - (Type(pTraits<Type>::one)*fluxStarAlbedo*gamma_/diffCoeff*3.0/4.0)/this->patch().deltaCoeffs()
	    );
    }

    fvPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::albedoSP3FvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const

{
    /*
    return tmp<Foam::Field<Type>>
    (
        new Foam::Field<Type>(this->size(), pTraits<Type>::zero)
    );
    */

    const Foam::Field<scalar>& diffCoeff =
        this->patch().template lookupPatchField<volScalarField, scalar>
        (
            diffCoeffName_
        );
    const Foam::Field<scalar>& fluxStarAlbedo =
        this->patch().template lookupPatchField<volScalarField, scalar>
        (
            fluxStarAlbedo_
        );

    if(forSecondMoment_)
    {
	    return Type(pTraits<Type>::one) * (1 + gamma_/diffCoeff*21.0/20.0/this->patch().deltaCoeffs());

    }
    else
    {
	    return Type(pTraits<Type>::one) * (1 + gamma_/diffCoeff/this->patch().deltaCoeffs());
    }


}

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::albedoSP3FvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const

{
    const Foam::Field<scalar>& diffCoeff =
        this->patch().template lookupPatchField<volScalarField, scalar>
        (
            diffCoeffName_
        );
    const Foam::Field<scalar>& fluxStarAlbedo =
        this->patch().template lookupPatchField<volScalarField, scalar>
        (
            fluxStarAlbedo_
        );

    if(forSecondMoment_)
    {
	    return -Type(pTraits<Type>::one) * (fluxStarAlbedo*gamma_/diffCoeff*3.0/20.0/this->patch().deltaCoeffs());

    }
    else
    {
	    return -Type(pTraits<Type>::one) * (fluxStarAlbedo*gamma_/diffCoeff*3.0/4.0/this->patch().deltaCoeffs());
    }


}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::albedoSP3FvPatchField<Type>::gradientInternalCoeffs() const
{
    const Foam::Field<scalar>& diffCoeff =
        this->patch().template lookupPatchField<volScalarField, scalar>
        (
            diffCoeffName_
        );
    const Foam::Field<scalar>& fluxStarAlbedo =
        this->patch().template lookupPatchField<volScalarField, scalar>
        (
            fluxStarAlbedo_
        );

    if(forSecondMoment_)
    {
	    return Type(pTraits<Type>::one)*(- gamma_/diffCoeff*21.0/20.0);

    }
    else
    {
	    return Type(pTraits<Type>::one)*(- gamma_/diffCoeff);
    }

}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::albedoSP3FvPatchField<Type>::gradientBoundaryCoeffs() const
{
    const Foam::Field<scalar>& diffCoeff =
        this->patch().template lookupPatchField<volScalarField, scalar>
        (
            diffCoeffName_
        );
    const Foam::Field<scalar>& fluxStarAlbedo =
        this->patch().template lookupPatchField<volScalarField, scalar>
        (
            fluxStarAlbedo_
        );

    if(forSecondMoment_)
    {
	    return Type(pTraits<Type>::one)*fluxStarAlbedo*gamma_/diffCoeff*3.0/20.0;

    }
    else
    {
	    return Type(pTraits<Type>::one)*fluxStarAlbedo*gamma_/diffCoeff*3.0/4.0;
    }
}

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::albedoSP3FvPatchField<Type>::snGrad() const
{
    const Foam::Field<scalar>& diffCoeff =
        this->patch().template lookupPatchField<volScalarField, scalar>
        (
            diffCoeffName_
        );
    const Foam::Field<scalar>& fluxStarAlbedo =
        this->patch().template lookupPatchField<volScalarField, scalar>
        (
            fluxStarAlbedo_
        );

    //return (- (this->patchInternalField()*gamma_/diffCoeff)/this->patch().deltaCoeffs());
    return (- (this->patchInternalField()*gamma_/diffCoeff));
}

template<class Type>
void Foam::albedoSP3FvPatchField<Type>::write(Ostream& os) const
/*{
    fvPatchField<Type>::write(os);
    gradient_.writeEntry("gradient", os);
}*/
{

    fvPatchField<Type>::write(os);

    os.writeKeyword("gamma")
        << gamma_ << token::END_STATEMENT << nl;
    os.writeKeyword("diffCoeffName")
        << diffCoeffName_ << token::END_STATEMENT << nl;
    os.writeKeyword("fluxStarAlbedo")
        << fluxStarAlbedo_ << token::END_STATEMENT << nl;
    os.writeKeyword("forSecondMoment")
        << forSecondMoment_ << token::END_STATEMENT << nl;

    os.writeKeyword("value uniform 0") << token::END_STATEMENT << nl;

    /*
    os.writeKeyword("gamma")
        << gamma_ << token::END_STATEMENT << nl;
    os.writeKeyword("diffCoeffName")
        << diffCoeffName_ << token::END_STATEMENT << nl;
    os.writeKeyword("fluxStarAlbedo")
        << fluxStarAlbedo_ << token::END_STATEMENT << nl;
    os.writeKeyword("forSecondMoment")
        << forSecondMoment_ << token::END_STATEMENT << nl;

    this->writeEntry("value", os);
    */

/*
    fvPatchField<Type>::write(os);
    //fixedGradientFvPatchScalarField::write(os);
    //writeEntryIfDifferent<word>(os, "diffCoeffName", "diffCoeffName", diffCoeffName_);


    gamma_.writeEntry("gamma", os);
    //diffCoeff_.writeEntry("diffCoeff", os);
    this->writeEntry("value", os);

    this->template writeEntryIfDifferent<word>(os, "diffCoeffName", "diffCoeffName", diffCoeffName_);
    this->template writeEntryIfDifferent<word>(os, "fluxStarAlbedo", "fluxStarAlbedo", fluxStarAlbedo_);
    this->template writeEntryIfDifferent<bool>(os, "forSecondMoment", "forSecondMoment", forSecondMoment_);
*/
}


// ************************************************************************* //
