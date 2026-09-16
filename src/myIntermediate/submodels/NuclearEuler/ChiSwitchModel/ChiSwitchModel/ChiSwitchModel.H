/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tommaso Pernatsch
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
    Foam::ChiSwitchModel

Group
    grpLagrangianIntermediateNuclearEulerSubModels

Description
    Base class for the rate constant (gamma) feeding
    NuclearEulerParcel::updateChi()'s per-Eulerian-step Bernoulli trial
    (chi switches with probability 1 - exp(-gamma*deltaT)).

    Unlike DispersionModel/ParticleForce, this is templated directly on
    NuclearEulerCloud<CloudType> (via CloudType::nuclearEulerCloudType),
    not on the fixed KinematicCloud layer - there is no pre-existing
    shared runtime-selection table to reuse here, so owner() can simply
    be given the full NuclearEulerCloud type instead, and every subclass
    gets direct, uncast access to alpha2()/U2()/d2()/k1()/epsilon1()/...

SourceFiles
    ChiSwitchModel.C

\*---------------------------------------------------------------------------*/

#ifndef ChiSwitchModel_H
#define ChiSwitchModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ChiSwitchModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ChiSwitchModel
:
    public CloudSubModelBase<CloudType>
{
public:

    //- Runtime type information
    TypeName("chiSwitchModel");


    // Declare runtime constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        ChiSwitchModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct from components
        ChiSwitchModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        ChiSwitchModel(const ChiSwitchModel<CloudType>& csm);

        //- Construct and return a clone
        virtual autoPtr<ChiSwitchModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~ChiSwitchModel();


    //- Selector
    static autoPtr<ChiSwitchModel<CloudType>> New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

        //- Return the switch rate constant gamma [1/s] for this parcel,
        //- at its current position (as recorded in trackingData)
        virtual scalar gamma
        (
            const typename CloudType::parcelType& p,
            const typename CloudType::parcelType::trackingData& td
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeChiSwitchModel(CloudType)                                        \
                                                                               \
    typedef Foam::CloudType::nuclearEulerCloudType nuclearEulerCloudType;    \
    defineTemplateTypeNameAndDebug                                           \
        (Foam::ChiSwitchModel<nuclearEulerCloudType>, 0);                    \
                                                                               \
    namespace Foam                                                           \
    {                                                                        \
        defineTemplateRunTimeSelectionTable                                  \
        (                                                                    \
            ChiSwitchModel<nuclearEulerCloudType>,                           \
            dictionary                                                       \
        );                                                                   \
    }

#define makeChiSwitchModelType(SS, CloudType)                                \
                                                                               \
    typedef Foam::CloudType::nuclearEulerCloudType nuclearEulerCloudType;    \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<nuclearEulerCloudType>, 0); \
                                                                               \
    Foam::ChiSwitchModel<nuclearEulerCloudType>::                            \
        adddictionaryConstructorToTable<Foam::SS<nuclearEulerCloudType>>     \
            add##SS##CloudType##nuclearEulerCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ChiSwitchModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
