/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::ForceMagnitudes

Group
    grpLagrangianIntermediateFunctionObjects

Description
    Diagnostic-only cloud function object. Writes the instantaneous
    magnitude of the drag, gravity+buoyancy, Saffman-Mei lift,
    pressure-gradient and virtual(added)-mass forces acting on every
    parcel, using the same formulas as the corresponding
    Foam::SphereDragForce, Foam::GravityForce, Foam::SaffmanMeiLiftForce,
    Foam::PressureGradientForce and Foam::VirtualMassForce classes.

    It does NOT feed back into particle motion (it only computes and
    writes force magnitudes after the real integration step has run), so
    adding it does not change parcel trajectories, regardless of which
    forces are active in the "particleForces" sub-dictionary.

    Operands:
    \table
      Operand        | Type         | Location
      output field   | scalarField  | \<time\>/lagrangian/\<cloud\>/FDrag
      output field   | scalarField  | \<time\>/lagrangian/\<cloud\>/FGrav
      output field   | scalarField  | \<time\>/lagrangian/\<cloud\>/FLift
      output field   | scalarField  | \<time\>/lagrangian/\<cloud\>/FPg
      output field   | scalarField  | \<time\>/lagrangian/\<cloud\>/FVm
    \endtable

Usage
    Minimal example by using \c constant/\<CloudProperties\>:
    \verbatim
    cloudFunctions
    {
        ForceMagnitudes1
        {
            // Mandatory entries
            type             ForceMagnitudes;
            rhoc             1758;
            nuc              1.5358e-06;

            // Optional entries
            U                U;
            Cvm              0.5;
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                              | Type   | Reqd | Deflt
      type      | Type name: ForceMagnitudes                | word   | yes  | -
      rhoc      | Carrier density                           | scalar | yes  | -
      nuc       | Carrier kinematic viscosity                | scalar | yes  | -
      U         | Name of carrier velocity field             | word   | no   | U
      Cvm       | Added-mass coefficient                     | scalar | no   | 0.5
    \endtable

    Note
      This model assumes a constant carrier density and viscosity
      (as supplied via rhoc/nuc), consistent with the Newtonian,
      constant-property transportModel used in this case. If the carrier
      properties vary in space/time in your case, replace the rhoc_/nuc_
      scalars with field lookups accordingly.

See also
  - Foam::SphereDragForce
  - Foam::GravityForce
  - Foam::SaffmanMeiLiftForce
  - Foam::PressureGradientForce
  - Foam::VirtualMassForce

SourceFiles
    ForceMagnitudes.C

\*---------------------------------------------------------------------------*/

#ifndef ForceMagnitudes_H
#define ForceMagnitudes_H

#include "CloudFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class ForceMagnitudes Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ForceMagnitudes
:
    public CloudFunctionObject<CloudType>
{
    // Private Data

        // Typedefs

            //- Convenience typedef for parcel type
            typedef typename CloudType::parcelType parcelType;

        //- Name of the carrier velocity field
        word UName_;

        //- Carrier density [kg/m3]
        scalar rhoc_;

        //- Carrier kinematic viscosity [m2/s]
        scalar nuc_;

        //- Added-mass coefficient
        scalar Cvm_;


public:

    //- Runtime type information
    TypeName("ForceMagnitudes");


    // Generated Methods

        //- No copy assignment
        void operator=(const ForceMagnitudes<CloudType>&) = delete;


    // Constructors

        //- Construct from dictionary
        ForceMagnitudes
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Copy construct
        ForceMagnitudes(const ForceMagnitudes<CloudType>& fm);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new ForceMagnitudes<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ForceMagnitudes() = default;


    // Member Functions

        //- Post-evolve hook
        virtual void postEvolve(const typename parcelType::trackingData& td);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ForceMagnitudes.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
