/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::DecayHeatModel

Group
    grpLagrangianIntermediateDecayHeatSubModels

Description
    Templated class to calculate the nuclear particle
    decay heat model.

Usage
    Minimal example by using \c constant/\<CloudProperties\>:
    \verbatim
    subModels
    {
        decayHeatModel    <model>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property           | Description                | Type | Reqd   | Deflt
      decayHeatModel     | Type name: \<model\>       | word | yes    | -
    \endtable

    Options for the \c \<model\> entry:
    \verbatim:
      constDecay      | Constant decay power model
      none            | No active model
    \endverbatim

SourceFiles
    DecayHeatModel.C
    DecayHeatModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef DecayHeatModel_H
#define DecayHeatModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DecayHeatModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class DecayHeatModel
:
    public CloudSubModelBase<CloudType>
{

public:

    //- Runtime type information
    TypeName("decayHeatModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        DecayHeatModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Generated Methods

        //- No copy assignment
        void operator=(const DecayHeatModel<CloudType>&) = delete;


    // Constructors

        //- Construct null from owner
        DecayHeatModel(CloudType& owner);

        //- Construct from dictionary
        DecayHeatModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Copy construct
        DecayHeatModel(const DecayHeatModel<CloudType>& dhm);

        //- Construct and return a clone
        virtual autoPtr<DecayHeatModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~DecayHeatModel() = default;


    //- Selector
    static autoPtr<DecayHeatModel<CloudType>> New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

        // Evaluation

            //- Return decay power based on particle diameter
            virtual scalar decayPower(const scalar d_part) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeDecayHeatModel(CloudType)                                       \
                                                                               \
    typedef Foam::CloudType::nuclearCloudType nuclearCloudType;                  \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::DecayHeatModel<nuclearCloudType>,                              \
        0                                                                      \
    );                                                                         \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            DecayHeatModel<nuclearCloudType>,                                \
            dictionary                                                         \
        );                                                                     \
    }


#define makeDecayHeatModelType(SS, CloudType)                               \
                                                                               \
    typedef Foam::CloudType::nuclearCloudType nuclearCloudType;                  \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<nuclearCloudType>, 0);         \
                                                                               \
    Foam::DecayHeatModel<nuclearCloudType>::                                 \
        adddictionaryConstructorToTable<Foam::SS<nuclearCloudType>>             \
            add##SS##CloudType##nuclearCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DecayHeatModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
