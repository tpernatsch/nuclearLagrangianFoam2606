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
    Foam::ConstantChiSwitchModel

Description
    Trivial chi-switch model: gamma is a single, dictionary-supplied
    constant for every parcel, everywhere, always - the original
    placeholder behaviour, kept as the default so existing cases (and
    chiSwitchGamma default 0, i.e. chi always false) are unaffected.

Usage
    \verbatim
    subModels
    {
        chiSwitchModel  constant;

        constantCoeffs
        {
            chiSwitchGamma  0;
        }
    }
    \endverbatim

SourceFiles
    ConstantChiSwitchModel.C

\*---------------------------------------------------------------------------*/

#ifndef ConstantChiSwitchModel_H
#define ConstantChiSwitchModel_H

#include "ChiSwitchModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class ConstantChiSwitchModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ConstantChiSwitchModel
:
    public ChiSwitchModel<CloudType>
{
    // Private Data

        //- The constant switch rate [1/s]
        const scalar gamma_;


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from components
        ConstantChiSwitchModel(const dictionary& dict, CloudType& owner);

        //- Construct copy
        ConstantChiSwitchModel(const ConstantChiSwitchModel<CloudType>& csm);

        //- Construct and return a clone
        virtual autoPtr<ChiSwitchModel<CloudType>> clone() const
        {
            return autoPtr<ChiSwitchModel<CloudType>>
            (
                new ConstantChiSwitchModel<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ConstantChiSwitchModel() = default;


    // Member Functions

        //- Return the switch rate constant gamma [1/s]
        virtual scalar gamma
        (
            const typename CloudType::parcelType& p,
            const typename CloudType::parcelType::trackingData& td
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ConstantChiSwitchModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
