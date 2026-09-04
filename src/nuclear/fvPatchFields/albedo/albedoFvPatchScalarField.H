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
    Foam::albedoFvPatchScalarField

Description
    This boundary condition applies an 'albedo' condition to neutron flux
    on an external boundary of the reactor, based on the 'mixed' type
    boundary condition:

    \f[
        \vec{n} \cdot \left( D \nabla \phi \right) + \frac{1}{2} \frac{1 - \beta}{1 + \beta} \phi = 0
    \f]

    where:
    \vartable
        \phi   | Neutron flux [m^-2 s^-1]
        D      | Diffusion coefficient [m]
        \beta  | Albedo coefficient [-]
    \endvartable

    The albedo coefficient \beta is defined as the ratio between the incoming
    and outcoming currents and can be set to any positive value. Particular
    cases are represented by 0 ('vacuum' or zero-incoming current condition)
    and 1 (symmetry condition).

Usage
    \table
    Property        | Description                      | Required   | Default value
    beta            | Albedo coefficient               | yes        |
    diffCoeffName   | Diffusion coefficient field name | no         | diffCoeffAlbedo
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            albedo;

        beta            uniform 0.0;

        value           $internalField;
    }
    \endverbatim

See also
    Foam::mixedFvPatchScalarField

SourceFiles
    albedoFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef albedoFvPatchScalarField_H
#define albedoFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class albedoFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class albedoFvPatchScalarField
:
    public mixedFvPatchScalarField
{
private:

    // Private Data

        //- Albedo coefficient
        scalarField beta_;

        //- Name of diffusion coefficient field
        const word diffCoeffName_;


public:

    //- Runtime type information
    TypeName("albedo");


    // Constructors

        //- Construct from patch and internal field
        albedoFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

       //- Construct from patch and internal field,
        //- initialise as zero-gradient
        albedoFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const Foam::zero
        );


	  //- Construct from patch, internal field and dictionary
        albedoFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  albedoFvPatchScalarField
        //  onto a new patch
        albedoFvPatchScalarField
        (
            const albedoFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        albedoFvPatchScalarField
        (
            const albedoFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new albedoFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        albedoFvPatchScalarField
        (
            const albedoFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new albedoFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Allow manipulation of the boundary values
            virtual bool fixesValue() const
            {
                return false;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap
		(
		const fvPatchFieldMapper&
		);


            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap
		(
		const fvPatchScalarField&,
		const labelList&
		);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
