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
    Foam::StickReplaceInteraction

Description
    Observes patch interactions and, for every parcel that an earlier
    patch interaction model in the chain has just stuck (kept in the
    cloud, marked inactive), clones it and injects a fresh replacement
    of the same size/properties at a random volume-weighted position in
    a named cell zone. Never deletes the stuck parcel, so the total
    parcel count grows by one per stick event - not a steady-state
    population (see "park and add" in the design discussion).

    Meant to be chained after the real stick/rebound/escape model via
    multiInteraction with oneInteractionOnly no - it only observes the
    keepParticle/active() state the earlier model already set, it never
    claims or alters the interaction itself.

Usage
    \verbatim
    patchInteractionModel multiInteraction;

    multiInteractionCoeffs
    {
        oneInteractionOnly no;

        model1
        {
            patchInteractionModel standardWallInteraction;
            standardWallInteractionCoeffs { type stick; ... }
        }
        model2
        {
            patchInteractionModel stickReplaceInteraction;
            stickReplaceInteractionCoeffs { cellZone liquid; }
        }
    }
    \endverbatim

SourceFiles
    StickReplaceInteraction.C

\*---------------------------------------------------------------------------*/

#ifndef StickReplaceInteraction_H
#define StickReplaceInteraction_H

#include "PatchInteractionModel.H"
#include "tetIndices.H"
#include "IDLList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class StickReplaceInteraction Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class StickReplaceInteraction
:
    public PatchInteractionModel<CloudType>
{
protected:

    typedef typename CloudType::parcelType parcelType;

    // Protected Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Name of the cell zone replacement parcels are placed in
        const word cellZoneName_;

        //- Cells of the zone on this processor
        labelList zoneCells_;

        //- Tet decomposition of each zone cell (for a volume-uniform
        //  random point inside it), and the cumulative tet-volume
        //  fraction within each cell used to pick which tet a point
        //  falls in - same technique as CellZoneInjection::setPositions
        List<List<tetIndices>> cellTetIs_;
        List<scalarList> cellTetCumulativeV_;

        //- Cumulative cell-volume fraction across this processor's own
        //  zone cells (0-1), used to pick which local cell a point
        //  falls in
        scalarList cellCumulativeV_;

        //- Cumulative zone volume per processor (raw, not normalised),
        //  size nProcs+1, index 0 = 0: processor i owns the range
        //  [sumProcV_[i], sumProcV_[i+1]) of the zone's total volume -
        //  same pattern as patchInjectionBase::sumTriMagSf_, but over
        //  cell volume instead of patch face area
        scalarList sumProcV_;

        //- Parcels stuck this step, cloned, waiting for postEvolve() to
        //  place their replacement
        IDLList<parcelType> stuckParcels_;

        //- Number/mass of replacement parcels injected so far
        label nInjected_;
        scalar massInjected_;


    // Protected Member Functions

        //- Cache the zone's cells, their tet decomposition and the
        //  per-processor cumulative volume table
        void buildZoneData();

        //- Processor owning the given fraction (0-1) of the zone's
        //  total volume, across all processors
        label whichProc(const scalar fraction01) const;

        //- A volume-uniform random position (and its owner cell/tet) in
        //  this processor's share of the zone
        void setPositionAndCell
        (
            Random& rnd,
            vector& position,
            label& cellOwner,
            label& tetFacei,
            label& tetPti
        ) const;

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);


public:

    //- Runtime type information
    TypeName("stickReplaceInteraction");


    // Constructors

        //- Construct from dictionary
        StickReplaceInteraction(const dictionary& dict, CloudType& cloud);

        //- Construct copy from owner cloud and patch interaction model
        StickReplaceInteraction(const StickReplaceInteraction<CloudType>& pim);

        //- Construct and return a clone using supplied owner cloud
        virtual autoPtr<PatchInteractionModel<CloudType>> clone() const
        {
            return autoPtr<PatchInteractionModel<CloudType>>
            (
                new StickReplaceInteraction<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~StickReplaceInteraction() = default;


    // Member Functions

        //- Observe the interaction: if the parcel has just been stuck
        //  (kept, but marked inactive) by an earlier model in the
        //  chain, clone it for replacement in postEvolve(). Never
        //  claims the interaction or changes keepParticle/active()
        //  itself.
        virtual bool correct
        (
            typename CloudType::parcelType& p,
            const polyPatch& pp,
            bool& keepParticle
        );

        //- Inject one replacement per parcel that stuck this step
        virtual void postEvolve();

        //- Write patch interaction info
        virtual void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "StickReplaceInteraction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
