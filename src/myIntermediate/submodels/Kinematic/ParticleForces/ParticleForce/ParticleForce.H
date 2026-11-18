/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::ParticleForce

Group
    grpLagrangianIntermediateForceSubModels

Description
    Abstract base class for particle forces

SourceFiles
    ParticleForce.C
    ParticleForceNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_ParticleForce_H
#define Foam_ParticleForce_H

#include "dictionary.H"
#include "forceSuSp.H"
#include "fvMesh.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ParticleForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticleForce
{
    // Private Data

        //- Reference to the owner cloud
        CloudType& owner_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Force coefficients dictionary
        const dictionary coeffs_;


public:

    //- Runtime type information
    TypeName("particleForce");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        ParticleForce,
        dictionary,
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        ),
        (owner, mesh, dict)
    );


    //- Convenience typedef for return type
    typedef VectorSpace<Vector<vector>, vector, 2> returnType;


    // Constructors

        //- Construct from mesh
        ParticleForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& forceType,
            const bool readCoeffs
        );

        //- Construct copy
        ParticleForce(const ParticleForce& pf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new ParticleForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ParticleForce();


    //- Selector
    static autoPtr<ParticleForce<CloudType>> New
    (
        CloudType& owner,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& forceType
    );


    // Member Functions

    // Access

        //- Return const access to the cloud owner
        const CloudType& owner() const noexcept { return owner_; }

        //- Return reference to the cloud owner
        CloudType& owner() noexcept { return owner_; }

        //- Return the mesh database
        const fvMesh& mesh() const noexcept { return mesh_; }

        //- Return the force coefficients dictionary
        const dictionary& coeffs() const noexcept { return coeffs_; }


        // Evaluation

            //- Cache fields
            virtual void cacheFields(const bool store);

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

            //- Calculate the non-coupled force
            virtual forceSuSp calcNonCoupled
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

            //- Return a direct random position displacement over this
            //- sub-step (e.g. overdamped Brownian motion), bypassing the
            //- velocity/force ODE entirely. Default: none. Only overridden
            //- by force models that opt into this treatment; see
            //- BrownianMotionForce's "displacement" switch.
            virtual vector calcDisplacement
            (
                const typename CloudType::parcelType& p,
                const typename CloudType::parcelType::trackingData& td,
                const scalar dt
            ) const
            {
                return Zero;
            }

            //- Return the diffusivity [m2/s] driving calcDisplacement(),
            //- for callers (e.g. the tracking sub-step clamp) that need to
            //- bound the step size before calcDisplacement() is actually
            //- evaluated. Default: zero (no displacement contribution).
            virtual scalar displacementD
            (
                const typename CloudType::parcelType& p,
                const typename CloudType::parcelType::trackingData& td
            ) const
            {
                return 0;
            }

            //- Whether this force model is in direct-displacement mode
            //- (i.e. calcDisplacement()/displacementD() are live). Cheap
            //- flag the tracking loop checks before paying for the extra
            //- setCellValues() a correct displacement draw needs at the
            //- start of a parcel's first sub-cycle - see KinematicParcel::
            //- move(). Default: false.
            virtual bool hasDisplacement() const
            {
                return false;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ParticleForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeParticleForceModel(CloudType)                                      \
                                                                               \
    typedef Foam::CloudType::kinematicCloudType kinematicCloudType;            \
    defineNamedTemplateTypeNameAndDebug                                        \
        (Foam::ParticleForce<kinematicCloudType>, 0);                          \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            ParticleForce<kinematicCloudType>,                                 \
            dictionary                                                         \
        );                                                                     \
    }


#define makeParticleForceModelType(SS, CloudType)                              \
                                                                               \
    typedef Foam::CloudType::kinematicCloudType kinematicCloudType;            \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<kinematicCloudType>, 0);      \
                                                                               \
    Foam::ParticleForce<kinematicCloudType>::                                  \
        adddictionaryConstructorToTable<Foam::SS<kinematicCloudType>>          \
            add##SS##CloudType##kinematicCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
