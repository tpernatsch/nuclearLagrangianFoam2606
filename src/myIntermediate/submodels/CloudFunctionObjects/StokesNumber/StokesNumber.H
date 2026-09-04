/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::StokesNumber

Group
    grpLagrangianIntermediateFunctionObjects

Description
    Calculates and writes particle Stokes number field on the cloud.

    \f[
        \mathrm{St}_p =
            \frac{\tau \, | u_0 |}{L_0}
    \f]

    \f[
        \tau =
            \frac{4\rho_p \, d_p^2}{3\mu_c C_D Re_p}
    
    \f]

    \vartable
        \mathrm{St}_p  | Particle Stokes number
        \tau           | Particle relaxation time
        u_c            | Velocity of the carrier
        L_0            | Characteristic length
        rho_p          | Density of the particle
        d_p            | Diameter of the particle
        mu_c           | Dynamic viscosity of carrier
    \endvartable

    Operands:
    \table
      Operand        | Type         | Location
      input          | -            | -
      output file    | -            | -
      output field   | scalarField  | \<time\>/lagrangian/\<cloud\>/Stk
    \endtable

Usage
    Minimal example by using \c constant/\<CloudProperties\>:
    \verbatim
    cloudFunctions
    {
        StokesNumber1
        {
            // Mandatory entries
            type             StokesNumber;
            L0               <value>;
            Uavg             <value>;
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                 | Type   | Reqd | Deflt
      type         | Type name: StokesNumber     | word   | yes  | -
      L0           | Type name: L0               | float  | yes  | -
      Uavg         | Type name: Uavg             | float  | yes  | -
    \endtable

See also
  - Foam::ThermoReynoldsNumber

SourceFiles
    StokesNumber.C

\*---------------------------------------------------------------------------*/

#ifndef StokesNumber_H
#define StokesNumber_H

#include "CloudFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class StokesNumber Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class StokesNumber
:
    public CloudFunctionObject<CloudType>
{
    // Private Data

        // Typedefs

            //- Convenience typedef for parcel type
            typedef typename CloudType::parcelType parcelType;


public:

    //- Runtime type information
    TypeName("StokesNumber");


    // Generated Methods

        scalar L0_;
        scalar Uavg_;
        
        //- No copy assignment
        void operator=(const StokesNumber<CloudType>&) = delete;


    // Constructors

        //- Construct from dictionary
        StokesNumber
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Copy construct
        StokesNumber(const StokesNumber<CloudType>& vf);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new StokesNumber<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~StokesNumber() = default;


    // Member Functions

        //- Post-evolve hook
        virtual void postEvolve(const typename parcelType::trackingData& td);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "StokesNumber.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
