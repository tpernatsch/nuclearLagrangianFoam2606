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
    Foam::NoDecayHeat

Group
    grpLagrangianIntermediateDecayHeatSubModels

Description
    Decay power model providing an interface to the
    properties of \c DecayHeatModel without any active model.

Usage
    Minimal example by using \c constant/\<CloudProperties\>:
    \verbatim
    subModels
    {
        decayHeatModel    none;
    }
    \endverbatim

    where the entries mean:
    \table
      Property           | Description              | Type | Reqd | Deflt
      decayHeatModel     | Type name: none          | word | yes  | -
    \endtable

SourceFiles
    NoDecayHeat.C

\*---------------------------------------------------------------------------*/

#ifndef NoDecayHeat_H
#define NoDecayHeat_H

#include "DecayHeatModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class NoDecayHeat Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoDecayHeat
:
    public DecayHeatModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");


    // Generated Methods

        //- No copy assignment
        void operator=(const NoDecayHeat&) = delete;


    // Constructors

        //- Construct from dictionary
        NoDecayHeat(const dictionary&, CloudType& owner);

        //- Copy construct
        NoDecayHeat(const NoDecayHeat<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<DecayHeatModel<CloudType>> clone() const
        {
            return autoPtr<DecayHeatModel<CloudType>>
            (
                new NoDecayHeat<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoDecayHeat() = default;


    // Member Functions

        //- Flag to indicate whether model activates heat transfer model
        virtual bool active() const
        {
            return false;
        }


        // Evaluation

            //- Return decay power per particle based on diameter
            virtual scalar decayPower(const scalar d_part) const
            {
                return 0.0;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NoDecayHeat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
