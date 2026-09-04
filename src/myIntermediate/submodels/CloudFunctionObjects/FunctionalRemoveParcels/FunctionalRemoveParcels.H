/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::FunctionalRemoveParcels

Group
    grpLagrangianIntermediateFunctionObjects

Description
    Removes parcels that hit user-specified face zone faces, with a
    probability given by a runtime-selectable Function1<scalar> of the
    parcel diameter, P(d).

    On every face crossing on a tracked zone, a uniform random number
    r in [0, 1) is drawn and the parcel is deleted if r < P(d); otherwise
    it survives the crossing unchanged. This lets the same face zone act
    as a partial/selective sink, e.g. modelling a screen or membrane whose
    capture efficiency depends on droplet size.

    Example usage:
    \verbatim
    functionalRemoveParcels1
    {
        type            functionalRemoveParcels;
        log             yes;
        resetOnWrite    no;
        resetOnStart    no;
        faceZones       (cycLeft cycRight);

        // Removal probability as a function of parcel diameter (m)
        probability
        {
            type        expression;
            expression  "0.5*(1 - tanh((arg() - dCrit)/dWidth))";
            variables
            (
                "dCrit  = 1e-4"
                "dWidth = 2e-5"
            );
        }

        // Or simply a constant probability:
        // probability
        // {
        //     type    constant;
        //     value   0.3;
        // }
    }
    \endverbatim

SourceFiles
    FunctionalRemoveParcels.C

\*---------------------------------------------------------------------------*/

#ifndef FunctionalRemoveParcels_H
#define FunctionalRemoveParcels_H

#include "CloudFunctionObject.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class OFstream;

/*---------------------------------------------------------------------------*\
                 Class FunctionalRemoveParcels Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class FunctionalRemoveParcels
:
    public CloudFunctionObject<CloudType>
{
        //- Convenience typedef for parcel type
        typedef typename CloudType::parcelType parcelType;


    // Private Data

        //- Face zone IDs
        labelList faceZoneIDs_;

        //- Number of parcels removed per face zone
        List<label> nParcels_;

        //- Mass of parcels removed per face zone
        List<scalar> mass_;

        //- Type ID of parcels to remove; all selected if not set (<0)
        label typeId_;

        //- Removal probability as a function of parcel diameter
        autoPtr<Function1<scalar>> probabilityPtr_;

        //- Flag to indicate whether data should be written to file
        bool logToFile_;

        //- Flag to reset counters on write
        bool resetOnWrite_;

        //- Flag to reset on start
        bool resetOnStart_;

        //- Output file pointer per zone
        PtrList<OFstream> outputFilePtr_;


    // Private Member Functions

        //- Helper function to create log files
        void makeLogFile
        (
            const word& zoneName,
            const label zoneI,
            const label nFaces,
            const scalar totArea
        );


protected:

    // Protected Member Functions

        //- Write post-processing info
        void write();


public:

    //- Runtime type information
    TypeName("functionalRemoveParcels");


    // Constructors

        //- Construct from dictionary
        FunctionalRemoveParcels
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        FunctionalRemoveParcels(const FunctionalRemoveParcels<CloudType>& ppm);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new FunctionalRemoveParcels<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~FunctionalRemoveParcels() = default;


    // Member Functions

        //- Post-evolve hook
        virtual void postEvolve
        (
            const typename parcelType::trackingData& td
        );

        //- Post-face hook
        virtual bool postFace
        (
            const parcelType& p,
            const typename parcelType::trackingData& td
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FunctionalRemoveParcels.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
