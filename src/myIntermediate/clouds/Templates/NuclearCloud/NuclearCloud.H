/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::NuclearCloud

Group
    grpLagrangianIntermediateClouds

Description
    Templated base class for thermodynamic nuclear cloud

    - Adds to kinematic cloud
      - Heat transfer
      - Decay heat
      - Decay population

SourceFiles
    NuclearCloudI.H
    NuclearCloud.C

\*---------------------------------------------------------------------------*/

#ifndef NuclearCloud_H
#define NuclearCloud_H

#include "KinematicCloud.H"
#include "nuclearCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class HeatTransferModel;

template<class CloudType>
class DecayHeatModel;

/*---------------------------------------------------------------------------*\
                        Class NuclearCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NuclearCloud
:
    public CloudType,
    public nuclearCloud
{
public:

    // Public typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::particleType parcelType;

        //- Convenience typedef for this cloud type
        typedef NuclearCloud<CloudType> nuclearCloudType;


private:

    // Private data

        //- Cloud copy pointer
        autoPtr<NuclearCloud<CloudType>> cloudCopyPtr_;


    // Private member functions

        //- No copy construct
        NuclearCloud(const NuclearCloud&) = delete;

        //- No copy assignment
        void operator=(const NuclearCloud&) = delete;


protected:

    // Protected data

        //- Thermo parcel constant properties
        typename parcelType::constantProperties constProps_;


        // References to the carrier fluid fields

            //- Vorticity [1/s]
            const volVectorField& vort_;

            //- Temperature [K]
            const volScalarField& T_;


        // References to the cloud sub-models

            //- Heat transfer model
            autoPtr<HeatTransferModel<NuclearCloud<CloudType>>>
                heatTransferModel_;

            //- Decay heat model
            autoPtr<DecayHeatModel<NuclearCloud<CloudType>>>
                decayHeatModel_;


        // Reference to the particle integration schemes

            //- Temperature integration
            autoPtr<integrationScheme> TIntegrator_;


        // Modelling options

            //- Include radiation
            Switch radiation_;

            //- Radiation sum of parcel projected areas
            autoPtr<volScalarField::Internal> radAreaP_;

            //- Radiation sum of parcel temperature^4
            autoPtr<volScalarField::Internal> radT4_;

            //- Radiation sum of parcel projected areas * temperature^4
            autoPtr<volScalarField::Internal> radAreaPT4_;


        // Sources

            //- Sensible enthalpy transfer [J/kg]
            autoPtr<volScalarField::Internal> hsTrans_;

            //- Coefficient for carrier phase temperature equation [J/K]
            autoPtr<volScalarField::Internal> hsCoeff_;

            //- Coefficient for carrier phase temperature equation multiplied by parcel temperature [J]
            autoPtr<volScalarField::Internal> hsCoeffTemp_;


    // Protected Member Functions

         // Initialisation

            //- Set cloud sub-models
            void setModels();


        // Cloud evolution functions

            //- Reset state of cloud
            void cloudReset(NuclearCloud<CloudType>& c);


public:

    // Constructors

        //- Construct given carrier gas fields
        NuclearCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& muc,
            const volVectorField& vort,
            const volScalarField& T,
            const dimensionedVector& g,
            bool readFields = true
        );

        //- Copy constructor with new name
        NuclearCloud(NuclearCloud<CloudType>& c, const word& name);

        //- Copy constructor with new name - creates bare cloud
        NuclearCloud
        (
            const fvMesh& mesh,
            const word& name,
            const NuclearCloud<CloudType>& c
        );

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> clone(const word& name)
        {
            return autoPtr<Cloud<parcelType>>
            (
                new NuclearCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType>>
            (
                new NuclearCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~NuclearCloud() = default;


    // Member Functions

        // Access

            //- Return a reference to the cloud copy
            inline const NuclearCloud& cloudCopy() const;

            //- Return the constant properties
            inline const typename parcelType::constantProperties&
                constProps() const;

            //- Return access to the constant properties
            inline typename parcelType::constantProperties& constProps();

            //- Return const access to the carrier vorticity field
            inline const volVectorField& vort() const;

            //- Return const access to the carrier temperature field
            inline const volScalarField& T() const;


            // Sub-models

                //- Return reference to heat transfer model
                inline const HeatTransferModel<NuclearCloud<CloudType>>&
                    heatTransfer() const;

                //- Return reference to decay heat model
                inline const DecayHeatModel<NuclearCloud<CloudType>>&
                    decayHeat() const;


            // Integration schemes

                //-Return reference to velocity integration
                inline const integrationScheme& TIntegrator() const;


            // Modelling options

                 //- Radiation flag
                inline bool radiation() const;

                //- Radiation sum of parcel projected areas [m2]
                inline volScalarField::Internal& radAreaP();

                //- Radiation sum of parcel projected areas [m2]
                inline const volScalarField::Internal&
                    radAreaP() const;

                //- Radiation sum of parcel temperature^4 [K4]
                inline volScalarField::Internal& radT4();

                //- Radiation sum of parcel temperature^4 [K4]
                inline const volScalarField::Internal& radT4() const;

                //- Radiation sum of parcel projected area*temperature^4 [m2K4]
                inline volScalarField::Internal& radAreaPT4();

                //- Radiation sum of parcel temperature^4 [m2K4]
                inline const volScalarField::Internal&
                    radAreaPT4() const;


            // Sources

                //- Transfer the effect of parcel to the carrier phase
                inline void transferToCarrier
                (
                    const parcelType& p,
                    const typename parcelType::trackingData& td
                );


                // Enthalpy

                    //- Sensible enthalpy transfer [J/kg]
                    inline volScalarField::Internal& hsTrans();

                    //- Sensible enthalpy transfer [J/kg]
                    inline const volScalarField::Internal&
                        hsTrans() const;

                    //- Return coefficient for carrier phase temperature equation
                    inline volScalarField::Internal& hsCoeff();

                    //- Return coefficient multiplied by parcel temperature for carrier phase temperature equation
                    inline volScalarField::Internal& hsCoeffTemp();

                    //- Return const coefficient for carrier phase temperature equation
                    inline const volScalarField::Internal&
                        hsCoeff() const;
                    
                    //- Return const coefficient multiplied by parcel temperature for carrier phase temperature equation
                    inline const volScalarField::Internal&
                        hsCoeffTemp() const;

                    //- Return sensible enthalpy source term [J/kg/m3/s]
                    inline tmp<fvScalarMatrix> Sh(volScalarField& hs, const dimensionedScalar& rho,const dimensionedScalar& Cp) const;


                // Radiation - overrides nuclearCloud virtual abstract members

                    //- Return tmp equivalent particulate emission
                    inline tmp<volScalarField> Ep() const;

                    //- Return tmp equivalent particulate absorption
                    inline tmp<volScalarField> ap() const;

                    //- Return tmp equivalent particulate scattering factor
                    inline tmp<volScalarField> sigmap() const;


        // Check

            //- Maximum temperature
            inline scalar Tmax() const;

            //- Minimum temperature
            inline scalar Tmin() const;


        // Cloud evolution functions

            //- Set parcel thermo properties
            void setParcelThermoProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt
            );

            //- Check parcel properties
            void checkParcelProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt,
                const bool fullyDescribed
            );

            //- Store the current cloud state
            void storeState();

            //- Reset the current cloud to the previously stored state
            void restoreState();

            //- Reset the cloud source terms
            void resetSourceTerms();

            //- Apply relaxation to (steady state) cloud sources
            void relaxSources(const NuclearCloud<CloudType>& cloudOldTime);

            //- Apply scaling to (transient) cloud sources
            void scaleSources();

            //- Pre-evolve
            void preEvolve(const typename parcelType::trackingData& td);

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

#include "NuclearCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NuclearCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
