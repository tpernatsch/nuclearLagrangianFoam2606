/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2019 OpenCFD Ltd.
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
    Foam::LocalInteractionMix

Group
    grpLagrangianIntermediatePatchInteractionSubModels

Description
    Patch interaction specified on a patch-by-patch basis

\*---------------------------------------------------------------------------*/

#ifndef LocalInteractionMix_H
#define LocalInteractionMix_H

#include "PatchInteractionModel.H"
#include "patchInteractionDataList.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class LocalInteractionMix Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class LocalInteractionMix
:
    public PatchInteractionModel<CloudType>
{
    // Private data

        // Reference to fluid mesh
        //const polyMesh& mesh_;

        //- List of participating patches
        const patchInteractionDataList patchData_;


        // Bookkeeping for particle fates

            //- Number of parcels escaped
            List<List<label>> nEscape_;

            //- Mass of parcels escaped
            List<List<scalar>> massEscape_;

            //- Number of parcels stuck to patches
            List<List<label>> nStick_;

            //- Mass of parcels stuck to patches
            List<List<scalar>> massStick_;


        //- Flag to output data as fields
        bool writeFields_;

        //- InjectorId to index map, when outputting escaped/stick/...
        //  particles sorted by injectorID
        Map<label> injIdToIndex_;

        //- Mass escape field
        autoPtr<volScalarField> massEscapePtr_;

        //- Mass stick field
        autoPtr<volScalarField> massStickPtr_;

        virtual void stick
        (
        typename CloudType::parcelType& p,
        const polyPatch& pp,
        bool& keepParticle,
        const label patchi,
        vector& U,
        const label idx
        );

        virtual void rebound
        (
        typename CloudType::parcelType& p,
        const polyPatch& pp,
        bool& keepParticle,
        const label patchi,
        vector& U,
        const label idx,
        const scalar e_n,
        const scalar e_t
        );

protected:

    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);


public:

    //- Runtime type information
    TypeName("localInteractionMix");


    // Constructors

        //- Construct from dictionary
        LocalInteractionMix(const dictionary& dict, CloudType& owner);

        //- Construct copy from owner cloud and patch interaction model
        LocalInteractionMix(const LocalInteractionMix<CloudType>& pim);

        //- Construct and return a clone using supplied owner cloud
        virtual autoPtr<PatchInteractionModel<CloudType>> clone() const
        {
            return autoPtr<PatchInteractionModel<CloudType>>
            (
                new LocalInteractionMix<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~LocalInteractionMix() = default;


    // Member Functions

        //- Return access to the massEscape field
        volScalarField& massEscape();

        //- Return access to the massStick field
        volScalarField& massStick();

        //- Apply velocity correction
        //  Returns true if particle remains in same cell
        virtual bool correct
        (
            typename CloudType::parcelType& p,
            const polyPatch& pp,
            bool& keepParticle
        );


        //- Write patch interaction info
        virtual void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LocalInteractionMix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
