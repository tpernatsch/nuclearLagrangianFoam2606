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
    Foam::depositionFvPatchScalarField

Description
    This boundary condition applies a 'deposition' condition to a scalar field,
    based on the 'mixed' type boundary condition:

    \f[
        - D \nabla c \cdot \vec{n} = \gamma (c - c_0)
    \f]

    where:
    \vartable
        c       | concentration [m^-3]
        D       | Diffusion coefficient [m^2 s^-1]
        \gamma  | deposition velocity [m s^-1]
        c_0     | reference concentration [m^-3]
    \endvartable

Usage
    \table
    Property        | Description                       | Required  | Default value
    gamma           | Deposition velocity               | yes       |
    c0              | Reference concentration           | yes       |
    diffCoeffName   | Diffusion coefficient field name  | no        | diffCoeffDeposition
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            deposition;

        gamma           uniform 0.0;
        c0              uniform 0.0;

        value           $internalField;
    }
    \endverbatim

See also
    Foam::mixedFvPatchScalarField

SourceFiles
    depositionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef depositionFvPatchScalarField_H
#define depositionFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class depositionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class depositionFvPatchScalarField
:
    public mixedFvPatchScalarField
{
private:

    // Private data

        //- Deposition velocity
        scalarField gamma_;

        //- Reference concentration
        scalarField c0_;

        //- Name of diffusion coefficient field
        const word diffCoeffName_;

public:

    //- Runtime type information
    TypeName("deposition");


    // Constructors

        //- Construct from patch and internal field
        depositionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

	//- Construct from patch and internal field,
        //- initialise as zero-gradient
        depositionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const Foam::zero
        );


        //- Construct from patch, internal field and dictionary
        depositionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  depositionFvPatchScalarField
        //  onto a new patch
        depositionFvPatchScalarField
        (
            const depositionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- //- Copy constructor
        depositionFvPatchScalarField
        (
            const depositionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new depositionFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        depositionFvPatchScalarField
        (
            const depositionFvPatchScalarField&,
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
                new depositionFvPatchScalarField(*this, iF)
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
            virtual void autoMap(const fvPatchFieldMapper&);


            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchScalarField&, const labelList&);



        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
