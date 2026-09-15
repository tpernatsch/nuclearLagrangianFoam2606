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
    Foam::NuclearEulerCloud

Group
    grpLagrangianIntermediateClouds

Description
    Templated Euler-Euler-aware nuclear cloud.

    Directly extends NuclearCloud (i.e. sits one layer above it in the
    Cloud template stack, exactly the way ReactingCloud sits above
    ThermoCloud, or ThermoCloud above KinematicCloud) so it inherits every
    NuclearCloud model unchanged: heat transfer, decay heat/decay
    population, radiation, and all of KinematicCloud's drag/dispersion/
    injection/patch-interaction machinery beneath that.

    What it adds is awareness of a *second* Euler-Euler carrier phase, on
    top of the single continuous-phase fields (rho, U, mu, T) that
    NuclearCloud/KinematicCloud already store. This is aimed at Lagrangian
    particle tracking inside a genuinely two-fluid (Euler-Euler) carrier,
    e.g. solid/nuclear particles tracked through the liquid phase of a
    gas-liquid bubble column, where the particle should still be able to
    "see" the gas phase's velocity, density, temperature and turbulence,
    and the local void fraction of each phase.

    Field roles:
      - The phase passed as (rho, U, muc, T) to the NuclearCloud
        base constructor is the "primary" or "tracking" phase: the
        particle physically moves through it, and all of the inherited
        drag/buoyancy/patch-interaction force models act against it,
        completely unchanged from plain NuclearCloud/KinematicCloud
        behaviour.
      - alpha1() is that primary phase's local volume fraction.
      - The bundle passed as "phase2" (see eulerCarrierPhaseFields)
        supplies the secondary phase's own volume fraction alpha2(),
        velocity U2(), density rho2(), viscosity mu2(), temperature T2()
        and turbulence closure fields k2()/epsilon2().
      - k1()/epsilon1() are the *primary* phase's own turbulence closure
        fields, stored here (rather than on KinematicCloud, which has no
        turbulence awareness at all) purely so that both phases can be
        combined consistently -- see NuclearEulerParcel.

    These extra fields are read (interpolated to the particle position)
    every sub-step by NuclearEulerParcel, which uses them for:
      - a void-fraction-weighted, two-phase convective heat transfer
        closure (the particle exchanges heat with both phases in
        proportion to their local volume fraction), and
      - a void-fraction-weighted turbulent-dispersion correction that
        adds extra stochastic agitation where the secondary phase is both
        locally present and more turbulent than the primary phase (e.g.
        bubble-induced turbulence acting on a particle nominally carried
        by the liquid).

    Momentum coupling (drag, buoyancy) is deliberately left against the
    primary phase only, using the existing, unmodified NuclearParcel/
    KinematicParcel force framework -- see the accompanying report for the
    reasoning. Two-way coupling back onto either Eulerian phase is out of
    scope (this cloud only *reads* the two-phase flow field); the usual
    one-way source terms inherited from NuclearCloud (hsTrans/hsCoeff/...)
    are still accumulated against the primary phase only, ready for a
    future extension.

SourceFiles
    NuclearEulerCloudI.H
    NuclearEulerCloud.C

\*---------------------------------------------------------------------------*/

#ifndef NuclearEulerCloud_H
#define NuclearEulerCloud_H

#include "NuclearCloud.H"
#include "nuclearEulerCloud.H"
#include "DispersionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class NuclearEulerCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NuclearEulerCloud
:
    public CloudType,
    public nuclearEulerCloud
{
public:

    // Public typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::particleType parcelType;

        //- Convenience typedef for this cloud type
        typedef NuclearEulerCloud<CloudType> nuclearEulerCloudType;


private:

    // Private data

        //- Cloud copy pointer
        autoPtr<NuclearEulerCloud<CloudType>> cloudCopyPtr_;


    // Private member functions

        //- No copy construct
        NuclearEulerCloud(const NuclearEulerCloud&) = delete;

        //- No copy assignment
        void operator=(const NuclearEulerCloud&) = delete;


protected:

    // Protected data

        // References to the secondary Euler carrier phase fields

            //- Primary (tracking) phase volume fraction [-]
            const volScalarField& alpha1_;

            //- Secondary phase volume fraction [-]
            const volScalarField& alpha2_;

            //- Secondary phase velocity [m/s]
            const volVectorField& U2_;

            //- Secondary phase density [kg/m3]
            const volScalarField& rho2_;

            //- Secondary phase dynamic viscosity [Pa.s]
            const volScalarField& mu2_;

            //- Secondary phase temperature [K]
            const volScalarField& T2_;

            //- Primary phase turbulence kinetic energy [m2/s2]
            const volScalarField& k1_;

            //- Primary phase turbulence dissipation rate [m2/s3]
            const volScalarField& epsilon1_;

            //- Secondary phase turbulence kinetic energy [m2/s2]
            const volScalarField& k2_;

            //- Secondary phase turbulence dissipation rate [m2/s3]
            const volScalarField& epsilon2_;

            //- True if the primary phase's turbulence model is RAS
            const bool k1IsRAS_;

            //- True if the secondary phase's turbulence model is RAS
            const bool k2IsRAS_;

            //- Dispersion model for the secondary phase
            autoPtr<DispersionModel<typename CloudType::kinematicCloudType>>
                dispersion2_;


public:

    // Constructors

        //- Construct given carrier fields for both Euler phases
        NuclearEulerCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& muc,
            const volScalarField& T,
            const dimensionedVector& g,
            const volScalarField& alpha1,
            const volScalarField& k1,
            const volScalarField& epsilon1,
            const bool k1IsRAS,
            const eulerCarrierPhaseFields& phase2,
            bool readFields = true
        );

        //- Copy constructor with new name
        NuclearEulerCloud(NuclearEulerCloud<CloudType>& c, const word& name);

        //- Copy constructor with new name - creates bare cloud
        NuclearEulerCloud
        (
            const fvMesh& mesh,
            const word& name,
            const NuclearEulerCloud<CloudType>& c
        );

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> clone(const word& name)
        {
            return autoPtr<Cloud<parcelType>>
            (
                new NuclearEulerCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType>>
            (
                new NuclearEulerCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~NuclearEulerCloud() = default;


    // Member Functions

        // Access

            //- Return a reference to the cloud copy
            inline const NuclearEulerCloud& cloudCopy() const;

            //- Return const access to the primary phase volume fraction
            inline const volScalarField& alpha1() const;

            //- Return const access to the secondary phase volume fraction
            inline const volScalarField& alpha2() const;

            //- Return const access to the secondary phase velocity
            inline const volVectorField& U2() const;

            //- Return const access to the secondary phase density
            inline const volScalarField& rho2() const;

            //- Return const access to the secondary phase viscosity
            inline const volScalarField& mu2() const;

            //- Return const access to the secondary phase temperature
            inline const volScalarField& T2() const;

            //- Return const access to the primary phase turbulence k
            inline const volScalarField& k1() const;

            //- Return const access to the primary phase turbulence epsilon
            inline const volScalarField& epsilon1() const;

            //- Return const access to the secondary phase turbulence k
            inline const volScalarField& k2() const;

            //- Return const access to the secondary phase turbulence
            //- epsilon
            inline const volScalarField& epsilon2() const;

            //- Return true if the primary phase's turbulence model is RAS
            inline bool k1IsRAS() const;

            //- Return true if the secondary phase's turbulence model is
            //- RAS
            inline bool k2IsRAS() const;

            //- Return const access to the secondary-phase dispersion
            //- model
            inline const DispersionModel<typename CloudType::kinematicCloudType>&
                dispersion2() const;

            //- Return non-const access to the secondary-phase dispersion
            //- model
            inline DispersionModel<typename CloudType::kinematicCloudType>&
                dispersion2();


        // Cloud evolution functions

            //- Pre-evolve
            void preEvolve(const typename parcelType::trackingData& td);

            //- Post-evolve
            void postEvolve(const typename parcelType::trackingData& td);

            //- Evolve the cloud
            void evolve();


        // Mapping

            //- Remap the cells of particles corresponding to the
            //  mesh topology change with a default tracking data object
            virtual void autoMap(const mapPolyMesh&);


        // I-O

            //- Print cloud information
            void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NuclearEulerCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NuclearEulerCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
