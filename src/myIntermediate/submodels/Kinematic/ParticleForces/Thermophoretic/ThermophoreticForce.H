/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::ThermophoreticForce

Description
    Calculates particle Thermophoretic force

SourceFiles
    ThermophoreticForceI.H
    ThermophoreticForce.C

\*---------------------------------------------------------------------------*/

#ifndef ThermophoreticForce_H
#define ThermophoreticForce_H

#include "ParticleForce.H"
#include "volFields.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class PressureGradientForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ThermophoreticForce
:
    public ParticleForce<CloudType>
{
protected:

    // Protected data

        //- Name of Temperature field
        const word   TName_;

        //- Particle thermal conductivity (W/mK) 
        const scalar Kp_;

        //- Carrier phase thermal conductivity (W/mK)
        const scalar Kc_;
  
        //- Carrier phase molecular mean free path (m)
        const scalar lambda_;

        //- Carrier phase temperature gradient interploator
        autoPtr<interpolation<vector>> gradTInterpPtr_;


public:

    //- Runtime type information
    TypeName("Thermophoretic");


    // Constructors

        //- Construct from mesh
        ThermophoreticForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& forceType = typeName
        );

        //- Construct copy
        ThermophoreticForce(const ThermophoreticForce& pgf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new ThermophoreticForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ThermophoreticForce();


    // Member Functions

        // Access

            inline const interpolation<vector>& gradTInterp() const;
    
            inline scalar Kp() const;

            inline scalar Kc() const;

            inline scalar lambda() const;

        // Evaluation

            //- Cache fields
            virtual void cacheFields(const bool store);

            //- Calculate the non-coupled force
            virtual forceSuSp calcCoupled
            (
                const typename CloudType::parcelType& p,
                const typename CloudType::parcelType::trackingData& td,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;       

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ThermophoreticForceI.H"

#ifdef NoRepository
    #include "ThermophoreticForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
