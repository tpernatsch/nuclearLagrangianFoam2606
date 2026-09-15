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
    Foam::NuclearEulerParcel

Group
    grpLagrangianIntermediateParcels

Description
    Nuclear parcel with awareness of a second Euler-Euler carrier phase.

    Directly extends NuclearParcel (which itself extends KinematicParcel),
    so drag, buoyancy, patch interaction, injection, decay heat and
    radiation are all inherited completely unchanged, evaluated against the
    "primary"/tracking phase exactly as for a plain NuclearParcel. This
    class adds no new persisted per-parcel state at all -- every extra
    field it reads (the secondary phase's alpha, U, rho, mu, T, k and
    epsilon, plus the primary phase's own alpha, k and epsilon) is
    transient, re-interpolated to the particle position every tracking
    sub-step via trackingData, exactly the way NuclearParcel/KinematicParcel
    already do for the primary phase's rho/U/mu/T.

    Convective heat transfer (calcHeatTransferTwoPhase()) is evaluated as
    a void-fraction-weighted mean of the same htc(Re,Pr) correlation
    applied once with the primary phase's slip Reynolds number/
    temperature and once with the secondary phase's. This reduces
    exactly to plain NuclearParcel::calcHeatTransfer() in the limit
    alpha2 -> 0.

    Turbulent dispersion is deliberately NOT touched here: the real
    dispersion models live in myTurbulence (DispersionRASModel and its
    derivatives), not myIntermediate, and this class relies entirely on
    whatever DispersionModel the cloud has configured via the inherited,
    unmodified ParcelType::calcDispersion() - myTurbulence would need its
    own two-phase-aware dispersion model to generalise that, which is
    out of scope here.

    A new per-parcel switch parameter, chi_ (a bool, always initialised to
    false), selects which phase's properties the *coupled* particle forces
    (drag, lift, gravity, virtual mass, pressure gradient, thermophoresis,
    Brownian motion) are evaluated against: F = F1*(1-chi) + F2*chi, where
    F1/F2 are the same force computed with the primary/secondary phase's
    fields respectively, and chi (implicitly converted to 0/1) always
    selects one or the other in full - the particle is either in the
    primary phase or the secondary phase, never a blend of both. The
    concrete two-phase force models themselves live in
    submodels/NuclearEuler/ParticleForces alongside a new
    makeNuclearEulerParcelForces registration macro - see the
    accompanying report.

SourceFiles
    NuclearEulerParcelI.H
    NuclearEulerParcel.C
    NuclearEulerParcelIO.C
    NuclearEulerParcelTrackingDataI.H

\*---------------------------------------------------------------------------*/

#ifndef NuclearEulerParcel_H
#define NuclearEulerParcel_H

#include "particle.H"
#include "demandDrivenEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class NuclearEulerParcel;

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const NuclearEulerParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                    Class NuclearEulerParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class NuclearEulerParcel
:
    public ParcelType
{
public:

    //- Class to hold two-phase-Euler nuclear particle constant properties
    class constantProperties
    :
        public ParcelType::constantProperties
    {
        // Private data

            //- Floor applied to (alpha1 + alpha2) as a weighting-average
            //- denominator
            demandDrivenEntry<scalar> alphaSmall_;

            //- Rate constant for the placeholder chi-switching model
            demandDrivenEntry<scalar> chiSwitchGamma_;


    public:

        // Constructors

            //- Null constructor
            constantProperties();

            //- Copy constructor
            constantProperties(const constantProperties& cp);

            //- Construct from dictionary
            constantProperties(const dictionary& parentDict);


        // Member functions

            //- Return const access to the void-fraction-sum floor
            inline scalar alphaSmall() const;

            //- Return const access to the chi-switch rate constant
            inline scalar chiSwitchGamma() const;
    };


    class trackingData
    :
        public ParcelType::trackingData
    {
    private:

        // Private data

            // Interpolators for the secondary Euler phase and the
            // primary phase's turbulence closure fields

                //- Primary phase volume fraction interpolator
                autoPtr<interpolation<scalar>> alpha1Interp_;

                //- Secondary phase volume fraction interpolator
                autoPtr<interpolation<scalar>> alpha2Interp_;

                //- Secondary phase velocity interpolator
                autoPtr<interpolation<vector>> U2Interp_;

                //- Secondary phase density interpolator
                autoPtr<interpolation<scalar>> rho2Interp_;

                //- Secondary phase viscosity interpolator
                autoPtr<interpolation<scalar>> mu2Interp_;

                //- Secondary phase temperature interpolator
                autoPtr<interpolation<scalar>> T2Interp_;

                //- Primary phase turbulence k interpolator
                autoPtr<interpolation<scalar>> k1Interp_;

                //- Primary phase turbulence epsilon interpolator
                autoPtr<interpolation<scalar>> epsilon1Interp_;

                //- Secondary phase turbulence k interpolator
                autoPtr<interpolation<scalar>> k2Interp_;

                //- Secondary phase turbulence epsilon interpolator
                autoPtr<interpolation<scalar>> epsilon2Interp_;


            // Cached values at the particle position

                //- Primary phase volume fraction [-]
                scalar alpha1c_;

                //- Secondary phase volume fraction [-]
                scalar alpha2c_;

                //- Secondary phase velocity [m/s]
                vector U2c_;

                //- Secondary phase density [kg/m3]
                scalar rho2c_;

                //- Secondary phase viscosity [Pa.s]
                scalar mu2c_;

                //- Secondary phase temperature [K]
                scalar T2c_;

                //- Primary phase turbulence kinetic energy [m2/s2]
                scalar k1c_;

                //- Primary phase turbulence dissipation rate [m2/s3]
                scalar epsilon1c_;

                //- Secondary phase turbulence kinetic energy [m2/s2]
                scalar k2c_;

                //- Secondary phase turbulence dissipation rate [m2/s3]
                scalar epsilon2c_;


    public:

        typedef typename ParcelType::trackingData::trackPart trackPart;

        // Constructors

            //- Construct from components
            template <class TrackCloudType>
            inline trackingData
            (
                const TrackCloudType& cloud,
                trackPart part = ParcelType::trackingData::tpLinearTrack
            );


        // Member functions

            //- Return const access to the primary phase volume fraction
            //- interpolator
            inline const interpolation<scalar>& alpha1Interp() const;

            //- Return const access to the secondary phase volume fraction
            //- interpolator
            inline const interpolation<scalar>& alpha2Interp() const;

            //- Return const access to the secondary phase velocity
            //- interpolator
            inline const interpolation<vector>& U2Interp() const;

            //- Return const access to the secondary phase density
            //- interpolator
            inline const interpolation<scalar>& rho2Interp() const;

            //- Return const access to the secondary phase viscosity
            //- interpolator
            inline const interpolation<scalar>& mu2Interp() const;

            //- Return const access to the secondary phase temperature
            //- interpolator
            inline const interpolation<scalar>& T2Interp() const;

            //- Return const access to the primary phase turbulence k
            //- interpolator
            inline const interpolation<scalar>& k1Interp() const;

            //- Return const access to the primary phase turbulence
            //- epsilon interpolator
            inline const interpolation<scalar>& epsilon1Interp() const;

            //- Return const access to the secondary phase turbulence k
            //- interpolator
            inline const interpolation<scalar>& k2Interp() const;

            //- Return const access to the secondary phase turbulence
            //- epsilon interpolator
            inline const interpolation<scalar>& epsilon2Interp() const;

            //- Return/access the primary phase volume fraction
            inline scalar alpha1c() const;
            inline scalar& alpha1c();

            //- Return/access the secondary phase volume fraction
            inline scalar alpha2c() const;
            inline scalar& alpha2c();

            //- Return/access the secondary phase velocity
            inline const vector& U2c() const;
            inline vector& U2c();

            //- Return/access the secondary phase density
            inline scalar rho2c() const;
            inline scalar& rho2c();

            //- Return/access the secondary phase viscosity
            inline scalar mu2c() const;
            inline scalar& mu2c();

            //- Return/access the secondary phase temperature
            inline scalar T2c() const;
            inline scalar& T2c();

            //- Return/access the primary phase turbulence k
            inline scalar k1c() const;
            inline scalar& k1c();

            //- Return/access the primary phase turbulence epsilon
            inline scalar epsilon1c() const;
            inline scalar& epsilon1c();

            //- Return/access the secondary phase turbulence k
            inline scalar k2c() const;
            inline scalar& k2c();

            //- Return/access the secondary phase turbulence epsilon
            inline scalar epsilon2c() const;
            inline scalar& epsilon2c();
    };


protected:

    // Protected data

        // Parcel properties

            //- Phase-switch parameter: false = primary phase, true =
            //- secondary phase, never a blend of both
            bool chi_;


    // Protected Member Functions

        //- Placeholder chi-switching model, for testing only
        template<class TrackCloudType>
        void updateChi(TrackCloudType& cloud);

        //- Calculate new particle temperature from a void-fraction
        //- weighted, two-phase convective heat transfer closure
        template<class TrackCloudType>
        scalar calcHeatTransferTwoPhase
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar dt,           // timestep
            const scalar Re1,          // Reynolds number vs. primary phase
            const scalar Re2,          // Reynolds number vs. secondary phase
            const scalar Pr,           // Prandtl number - surface
            const scalar kappa,        // Thermal conductivity - surface
            const scalar NCpW,         // Sum of N*Cp*W of emission species
            const scalar Sh,           // explicit particle enthalpy source
            scalar& dhsTrans,          // sensible enthalpy transfer to carrier
            scalar& Sph                // linearised heat transfer coefficient
        );


public:

    //- Size in bytes of the fields
    static const std::size_t sizeofFields;


    // Static data members

        //- Runtime type information
        TypeName("NuclearEulerParcel");

        //- String representation of properties
        AddToPropertyList(ParcelType, " chi");


    // Constructors

        //- Construct from mesh, coordinates and topology
        //  Other properties initialised as null
        inline NuclearEulerParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti
        );

        //- Construct from a position and a cell, searching for the rest of
        //  the required topology. Other properties are initialised as null.
        inline NuclearEulerParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli
        );

        //- Construct from components
        inline NuclearEulerParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const scalar dTarget0,
            const vector& U0,
            const vector& f0,
            const vector& angularMomentum0,
            const vector& torque0,
            const constantProperties& constProps
        );

        //- Construct from Istream
        NuclearEulerParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct as a copy
        NuclearEulerParcel(const NuclearEulerParcel& p);

        //- Construct as a copy
        NuclearEulerParcel(const NuclearEulerParcel& p, const polyMesh& mesh);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new NuclearEulerParcel(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new NuclearEulerParcel(*this, mesh));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<NuclearEulerParcel<ParcelType>> operator()(Istream& is)
            const
            {
                return autoPtr<NuclearEulerParcel<ParcelType>>
                (
                    new NuclearEulerParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to the phase-switch parameter chi
            inline bool chi() const;


        // Edit

            //- Return access to the phase-switch parameter chi
            inline bool& chi();


    // Main calculation loop

        //- Set cell values
        template<class TrackCloudType>
        void setCellValues(TrackCloudType& cloud, trackingData& td);

        //- Apply dispersion, primary phase always plus secondary phase
        //- when switched and RAS
        template<class TrackCloudType>
        void calcDispersion
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar dt
        );

        //- Update parcel properties over the time interval
        template<class TrackCloudType>
        void calc
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar dt
        );


    // I-O

        //- Read
        template<class CloudType>
        static void readFields(CloudType& c);

        //- Write
        template<class CloudType>
        static void writeFields(const CloudType& c);

        //- Write individual parcel properties to stream
        void writeProperties
        (
            Ostream& os,
            const wordRes& filters,
            const word& delim,
            const bool namesOnly = false
        ) const;

        //- Read particle fields as objects from the obr registry
        template<class CloudType>
        static void readObjects(CloudType& c, const objectRegistry& obr);

        //- Write particle fields as objects into the obr registry
        template<class CloudType>
        static void writeObjects(const CloudType& c, objectRegistry& obr);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const NuclearEulerParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NuclearEulerParcelI.H"
#include "NuclearEulerParcelTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NuclearEulerParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
