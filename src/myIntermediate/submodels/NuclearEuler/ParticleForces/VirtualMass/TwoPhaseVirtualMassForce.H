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
    Foam::TwoPhaseVirtualMassForce

Group
    grpLagrangianIntermediateForceSubModels

Description
    Two-phase-Euler-aware particle virtual (added) mass force.

    Only usable with a cloud whose parcelType::trackingData exposes the
    secondary-phase/chi accessors added by NuclearEulerParcel.

    Scales the already phase-1/phase-2-blended
    TwoPhasePressureGradientForce result by a single virtual-mass
    coefficient Cvm (typically 0.5), and adds a matching chi-blended
    contribution to the effective (added) mass used by the velocity
    integrator:
    \f[
        \mathrm{massAdd} =
            m\left(\frac{\rho_{1c}}{\rho_p}(1-\chi)
                 + \frac{\rho_{2c}}{\rho_p}\chi\right) C_\mathrm{vm}
    \f]
    chi is currently always 0 (no switching model exists yet), so this
    reduces exactly to the primary-phase-only VirtualMassForce result.

Usage
    \verbatim
    subModels
    {
        particleForces
        {
            twoPhaseVirtualMass;
        }
    }

    twoPhaseVirtualMassCoeffs
    {
        Cvm 0.5;
        U1  U.water;   // optional, defaults to "U.water"
        U2  U.air;     // optional, defaults to "U.air"
    }
    \endverbatim

SourceFiles
    TwoPhaseVirtualMassForce.C

\*---------------------------------------------------------------------------*/

#ifndef TwoPhaseVirtualMassForce_H
#define TwoPhaseVirtualMassForce_H

#include "TwoPhasePressureGradientForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class TwoPhaseVirtualMassForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class TwoPhaseVirtualMassForce
:
    public TwoPhasePressureGradientForce<CloudType>
{
    // Private data

        //- Virtual mass coefficient - typically 0.5
        scalar Cvm_;


public:

    //- Runtime type information
    TypeName("twoPhaseVirtualMass");


    // Constructors

        //- Construct from mesh
        TwoPhaseVirtualMassForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& forceType = typeName
        );

        //- Construct copy
        TwoPhaseVirtualMassForce(const TwoPhaseVirtualMassForce& vmf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new TwoPhaseVirtualMassForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~TwoPhaseVirtualMassForce() = default;


    // Member Functions

        // Evaluation

            //- Calculate the coupled force
            virtual forceSuSp calcCoupled
            (
                const typename CloudType::parcelType& p,
                const typename CloudType::parcelType::trackingData& td,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;

            //- Return the added mass
            virtual scalar massAdd
            (
                const typename CloudType::parcelType& p,
                const typename CloudType::parcelType::trackingData& td,
                const scalar mass
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TwoPhaseVirtualMassForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
