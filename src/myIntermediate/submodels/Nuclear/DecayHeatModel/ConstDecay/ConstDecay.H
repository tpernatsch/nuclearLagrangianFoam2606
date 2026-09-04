/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::ConstDecay

Group
    grpLagrangianIntermediateDecayHeatSubModels

Description
    Constant decay heat model for nuclear particles.
    Decay power is calculated based on the particle
    diameter using a reference diameter and reference
    decay power. Moreover, a number of particles per
    parcel is assumed to be 100, to understand the 
    percentage of decay in each parcel.


Usage
    Minimal example by using \c constant/\<CloudProperties\>:
    \verbatim
    subModels
    {
        // Mandatory entries
        DecayHeatModel    ConstDecay;

    }
    \endverbatim

    where the entries mean:
    \table
      Property           | Description              | Type | Reqd | Deflt
      DecayHeatModel  | Type name: ConstDecay       | word | yes  | -
    \endtable

SourceFiles
    ConstDecay.C

\*---------------------------------------------------------------------------*/

#ifndef ConstDecay_H
#define ConstDecay_H

#include "DecayHeatModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ConstDecay Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ConstDecay
:
    public DecayHeatModel<CloudType>
{
    //private data
    
    scalar decayPower_;

    // Reference diameter for noble metal fission product [m]
    scalar d0_;


    // Reference decay power for noble metal fission product [W/particle]
    scalar Q0_;

public:

    //- Runtime type information
    TypeName("constDecay");


    // Generated Methods

        //- No copy assignment
        void operator=(const ConstDecay&) = delete;


    // Constructors

        //- Construct from dictionary
        ConstDecay(const dictionary& dict, CloudType& cloud);

        //- Copy construct
        ConstDecay(const ConstDecay<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<DecayHeatModel<CloudType>> clone() const
        {
            return autoPtr<DecayHeatModel<CloudType>>
            (
                new ConstDecay<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ConstDecay() = default;


    // Member Functions

        // Evaluation

            //- Return decay power based on particle diameter
            virtual scalar decayPower(const scalar d_part) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ConstDecay.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
