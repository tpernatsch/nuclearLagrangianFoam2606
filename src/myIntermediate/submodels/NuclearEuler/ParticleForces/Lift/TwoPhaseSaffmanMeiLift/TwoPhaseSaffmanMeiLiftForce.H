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
    Foam::TwoPhaseSaffmanMeiLiftForce

Group
    grpLagrangianIntermediateForceSubModels

Description
    Two-phase-Euler-aware Saffman-Mei lift force - see TwoPhaseLiftForce
    for the phase-blending mechanics and SaffmanMeiLiftForce for the
    correlation itself, applied here once per Euler phase.

Usage
    \verbatim
    subModels
    {
        particleForces
        {
            twoPhaseSaffmanMeiLift;
        }
    }
    \endverbatim

SourceFiles
    TwoPhaseSaffmanMeiLiftForce.C

\*---------------------------------------------------------------------------*/

#ifndef TwoPhaseSaffmanMeiLiftForce_H
#define TwoPhaseSaffmanMeiLiftForce_H

#include "TwoPhaseLiftForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class TwoPhaseSaffmanMeiLiftForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class TwoPhaseSaffmanMeiLiftForce
:
    public TwoPhaseLiftForce<CloudType>
{
protected:

    // Protected Member Functions

        //- Calculate the lift coefficient
        virtual scalar Cl
        (
            const typename CloudType::parcelType& p,
            const typename CloudType::parcelType::trackingData& td,
            const vector& curlUc,
            const scalar rhoc,
            const scalar Re,
            const scalar muc
        ) const;


public:

    //- Runtime type information
    TypeName("twoPhaseSaffmanMeiLift");


    // Constructors

        //- Construct from mesh
        TwoPhaseSaffmanMeiLiftForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& forceType = typeName
        );

        //- Construct copy
        TwoPhaseSaffmanMeiLiftForce(const TwoPhaseSaffmanMeiLiftForce& lf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new TwoPhaseSaffmanMeiLiftForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~TwoPhaseSaffmanMeiLiftForce() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TwoPhaseSaffmanMeiLiftForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
