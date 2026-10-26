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

\*---------------------------------------------------------------------------*/

#include "StickReplaceInteraction.H"
#include "polyMeshTetDecomposition.H"
#include "PstreamBuffers.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::StickReplaceInteraction<CloudType>::buildZoneData()
{
    const label zoneI = mesh_.cellZones().findZoneID(cellZoneName_);

    if (zoneI < 0)
    {
        FatalErrorInFunction
            << "Unknown cell zone name: " << cellZoneName_
            << ". Valid cell zones are: " << mesh_.cellZones().names()
            << nl << exit(FatalError);
    }

    zoneCells_ = mesh_.cellZones()[zoneI];

    const scalarField& V = mesh_.V();

    cellTetIs_.setSize(zoneCells_.size());
    cellTetCumulativeV_.setSize(zoneCells_.size());

    scalarList cellV(zoneCells_.size());

    forAll(zoneCells_, i)
    {
        const label celli = zoneCells_[i];
        cellV[i] = V[celli];

        cellTetIs_[i] = polyMeshTetDecomposition::cellTetIndices(mesh_, celli);

        const List<tetIndices>& tets = cellTetIs_[i];
        scalarList& vf = cellTetCumulativeV_[i];
        vf.setSize(tets.size(), Zero);

        for (label t = 1; t < tets.size() - 1; ++t)
        {
            vf[t] = vf[t-1] + tets[t].tet(mesh_).mag()/cellV[i];
        }
        if (tets.size())
        {
            vf.last() = 1.0;
        }
    }

    // Cumulative cell-volume fraction across this processor's own zone
    // cells, used to pick which local cell a point falls in
    const scalar localV = sum(cellV);

    cellCumulativeV_.setSize(zoneCells_.size(), Zero);
    scalar run = 0;
    forAll(zoneCells_, i)
    {
        run += cellV[i];
        cellCumulativeV_[i] = (localV > 0 ? run/localV : 0);
    }

    // Cumulative zone volume per processor (patchInjectionBase pattern)
    const label myProci = UPstream::myProcNo();
    const label numProc = UPstream::nProcs();

    sumProcV_.resize_nocopy(numProc + 1);
    sumProcV_[0] = 0;
    {
        scalarList::subList slice(sumProcV_, numProc, 1);
        slice[myProci] = localV;
        Pstream::allGatherList(slice);
    }

    for (label i = 1; i < sumProcV_.size(); ++i)
    {
        sumProcV_[i] += sumProcV_[i-1];
    }
}


template<class CloudType>
Foam::label Foam::StickReplaceInteraction<CloudType>::whichProc
(
    const scalar fraction01
) const
{
    const scalar volFraction = fraction01*sumProcV_.last();

    forAllReverse(sumProcV_, i)
    {
        if (volFraction >= sumProcV_[i])
        {
            return i;
        }
    }

    return 0;
}


template<class CloudType>
void Foam::StickReplaceInteraction<CloudType>::setPositionAndCell
(
    Random& rnd,
    vector& position,
    label& cellOwner,
    label& tetFacei,
    label& tetPti
) const
{
    // Pick a zone cell on this processor, weighted by volume
    const scalar cellFraction = rnd.sample01<scalar>();

    label i = cellCumulativeV_.size() - 1;
    forAll(cellCumulativeV_, ci)
    {
        if (cellCumulativeV_[ci] > cellFraction)
        {
            i = ci;
            break;
        }
    }

    // Pick a tet within that cell, weighted by volume, then a random
    // point within the tet
    const List<tetIndices>& tets = cellTetIs_[i];
    const scalarList& vf = cellTetCumulativeV_[i];

    const scalar tetFraction = rnd.sample01<scalar>();
    label tetI = tets.size() - 1;
    forAll(vf, t)
    {
        if (vf[t] > tetFraction)
        {
            tetI = t;
            break;
        }
    }

    position = tets[tetI].tet(mesh_).randomPoint(rnd);
    cellOwner = zoneCells_[i];
    tetFacei = tets[tetI].face();
    tetPti = tets[tetI].tetPt();
}


template<class CloudType>
void Foam::StickReplaceInteraction<CloudType>::writeFileHeader(Ostream& os)
{
    PatchInteractionModel<CloudType>::writeFileHeader(os);

    this->writeTabbed(os, cellZoneName_ + "_nInjected");
    this->writeTabbed(os, cellZoneName_ + "_massInjected");
}


template<class CloudType>
void Foam::StickReplaceInteraction<CloudType>::writeHits()
{
    // Collect the hits of all processors on the master (collective call)
    const List<scalar> t(globalIndex::gatherOp(hitTime_));
    const List<point> x(globalIndex::gatherOp(hitPosition_));
    const List<scalar> age(globalIndex::gatherOp(hitAge_));
    const labelList patchi(globalIndex::gatherOp(hitPatch_));
    const labelList typeId(globalIndex::gatherOp(hitTypeId_));
    const labelList origProc(globalIndex::gatherOp(hitOrigProc_));
    const labelList origId(globalIndex::gatherOp(hitOrigId_));

    hitTime_.clear();
    hitPosition_.clear();
    hitAge_.clear();
    hitPatch_.clear();
    hitTypeId_.clear();
    hitOrigProc_.clear();
    hitOrigId_.clear();

    if (!Pstream::master())
    {
        return;
    }

    forAll(t, i)
    {
        const word key
        (
            mesh_.boundaryMesh()[patchi[i]].name() + "_" + Foam::name(typeId[i])
        );

        label fi = hitFileIndex_.lookup(key, -1);

        if (fi < 0)
        {
            // One file per patch and injector, in the start-time folder
            autoPtr<OFstream> osPtr = this->newFileAtStartTime("stickHits_" + key);

            if (!osPtr)
            {
                continue;
            }

            osPtr->precision(10);
            *osPtr
                << "# Parcels stuck on " << mesh_.boundaryMesh()[patchi[i]].name()
                << ", injector (typeId) " << typeId[i] << nl
                << "# Time\tx\ty\tz\tage\torigProc\torigId" << nl;

            fi = hitFiles_.size();
            hitFiles_.resize(fi + 1);
            hitFiles_.set(fi, osPtr.release());
            hitFileIndex_.insert(key, fi);
        }

        hitFiles_[fi]
            << t[i] << tab
            << x[i].x() << tab << x[i].y() << tab << x[i].z() << tab
            << age[i] << tab << origProc[i] << tab << origId[i] << nl;
    }

    forAll(hitFiles_, fi)
    {
        hitFiles_[fi].flush();
    }

    nHitsWritten_ += t.size();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::StickReplaceInteraction<CloudType>::StickReplaceInteraction
(
    const dictionary& dict,
    CloudType& cloud
)
:
    PatchInteractionModel<CloudType>(dict, cloud, typeName),
    mesh_(cloud.mesh()),
    cellZoneName_(this->coeffDict().template get<word>("cellZone")),
    deleteStuck_(this->coeffDict().getOrDefault("deleteStuck", false)),
    writeHits_(this->coeffDict().getOrDefault("writeHits", false)),
    nInjected_(0),
    massInjected_(0),
    nHitsWritten_(0)
{
    buildZoneData();
}


template<class CloudType>
Foam::StickReplaceInteraction<CloudType>::StickReplaceInteraction
(
    const StickReplaceInteraction<CloudType>& pim
)
:
    PatchInteractionModel<CloudType>(pim),
    mesh_(pim.mesh_),
    cellZoneName_(pim.cellZoneName_),
    zoneCells_(pim.zoneCells_),
    cellTetIs_(pim.cellTetIs_),
    cellTetCumulativeV_(pim.cellTetCumulativeV_),
    cellCumulativeV_(pim.cellCumulativeV_),
    sumProcV_(pim.sumProcV_),
    deleteStuck_(pim.deleteStuck_),
    writeHits_(pim.writeHits_),
    hitTime_(pim.hitTime_),
    hitPosition_(pim.hitPosition_),
    hitAge_(pim.hitAge_),
    hitPatch_(pim.hitPatch_),
    hitTypeId_(pim.hitTypeId_),
    hitOrigProc_(pim.hitOrigProc_),
    hitOrigId_(pim.hitOrigId_),
    hitFiles_(),        // output files are not copied
    hitFileIndex_(),
    nInjected_(pim.nInjected_),
    massInjected_(pim.massInjected_),
    nHitsWritten_(pim.nHitsWritten_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::StickReplaceInteraction<CloudType>::correct
(
    typename CloudType::parcelType& p,
    const polyPatch& pp,
    bool& keepParticle
)
{
    // Only act on a parcel an earlier model in the chain has just parked
    // (kept, but marked inactive) - i.e. a genuine "stick" outcome, not
    // escape (keepParticle false) or rebound (still active)
    if (keepParticle && !p.active())
    {
        stuckParcels_.append(static_cast<parcelType*>(p.clone().ptr()));

        if (writeHits_)
        {
            // Time of the hit within the step, position on the wall face
            const Time& runTime = this->owner().db().time();

            hitTime_.push_back
            (
                runTime.value() - (1 - p.stepFraction())*runTime.deltaTValue()
            );
            hitPosition_.push_back(p.position());
            hitAge_.push_back(p.age());
            hitPatch_.push_back(pp.index());
            hitTypeId_.push_back(p.typeId());
            hitOrigProc_.push_back(p.origProc());
            hitOrigId_.push_back(p.origId());
        }

        if (deleteStuck_)
        {
            // Counted by the stick model and recorded above: remove it, so
            // the cloud does not grow with the deposits
            keepParticle = false;
        }
    }

    // Never claims the interaction (the stick model did)
    return false;
}


template<class CloudType>
void Foam::StickReplaceInteraction<CloudType>::postEvolve()
{
    auto& rnd = this->owner().rndGen();

    if (Pstream::parRun())
    {
        // See comments in Cloud::move() about transfer particles handling
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
        PtrList<UOPstream> UOPstreamPtrs(Pstream::nProcs());

        for (parcelType& p : stuckParcels_)
        {
            // Choose a random location to place the replacement, and the
            // processor that owns it
            const scalar fraction01 = rnd.template sample01<scalar>();
            const label toProci = whichProc(fraction01);

            auto* osptr = UOPstreamPtrs.get(toProci);
            if (!osptr)
            {
                osptr = new UOPstream(toProci, pBufs);
                UOPstreamPtrs.set(toProci, osptr);
            }

            (*osptr) << p;

            delete(stuckParcels_.remove(&p));
        }

        pBufs.finishedSends();

        for (const int proci : pBufs.allProcs())
        {
            if (pBufs.recvDataCount(proci))
            {
                UIPstream is(proci, pBufs);

                while (!is.eof())
                {
                    auto* newp = new parcelType(this->owner().mesh(), is);

                    vector newPosition;
                    label cellOwner, tetFacei, tetPti;
                    setPositionAndCell
                    (
                        rnd, newPosition, cellOwner, tetFacei, tetPti
                    );

                    newp->relocate(newPosition, cellOwner);
                    newp->active(true);
                    newp->age() = 0;
                    newp->U() = this->owner().U()[cellOwner];

                    // A new parcel, not the stuck one: own identity
                    // (typeId, the size group, is kept)
                    newp->origProc() = Pstream::myProcNo();
                    newp->origId() = newp->getNewParticleID();

                    ++nInjected_;
                    massInjected_ += newp->nParticle()*newp->mass();

                    this->owner().addParticle(newp);
                }
            }
        }
    }
    else
    {
        for (parcelType& p : stuckParcels_)
        {
            parcelType* newp = stuckParcels_.remove(&p);

            vector newPosition;
            label cellOwner, tetFacei, tetPti;
            setPositionAndCell(rnd, newPosition, cellOwner, tetFacei, tetPti);

            newp->relocate(newPosition, cellOwner);
            newp->active(true);
            newp->age() = 0;
            newp->U() = this->owner().U()[cellOwner];

            // A new parcel, not the stuck one: own identity
            // (typeId, the size group, is kept)
            newp->origProc() = Pstream::myProcNo();
            newp->origId() = newp->getNewParticleID();

            ++nInjected_;
            massInjected_ += newp->nParticle()*newp->mass();

            this->owner().addParticle(newp);
        }
    }

    // Hits to file at the write times (consistent with a restart from them)
    if (writeHits_ && this->owner().db().time().writeTime())
    {
        writeHits();
    }
}


template<class CloudType>
void Foam::StickReplaceInteraction<CloudType>::info()
{
    PatchInteractionModel<CloudType>::info();

    label npr0 = this->template getModelProperty<label>("nInjected");
    scalar mpr0 = this->template getModelProperty<scalar>("massInjected");

    label npr = nInjected_;
    reduce(npr, sumOp<label>());
    npr += npr0;

    scalar mpr = massInjected_;
    reduce(mpr, sumOp<scalar>());
    mpr += mpr0;

    Log_<< "    Stuck-parcel replacement: injected into " << cellZoneName_
        << " = " << npr << ", " << mpr << nl;

    if (writeHits_)
    {
        Log_<< "    Stuck-parcel hits written to file (this run) = "
            << returnReduce(nHitsWritten_, maxOp<label>()) << nl;
    }

    this->file() << tab << npr << tab << mpr << endl;

    if (this->writeTime())
    {
        this->setModelProperty("nInjected", npr);
        this->setModelProperty("massInjected", mpr);
        nInjected_ = 0;
        massInjected_ = 0;
    }
}


// ************************************************************************* //
