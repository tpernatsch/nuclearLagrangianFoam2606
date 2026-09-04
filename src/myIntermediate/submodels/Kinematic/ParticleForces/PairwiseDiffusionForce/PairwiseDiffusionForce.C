/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "PairwiseDiffusionForce.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolationCellPoint.H"
#include "Cloud.H"
#include "PstreamBuffers.H"

namespace Foam
{

template<class CloudType>
PairwiseDiffusionForce<CloudType>::PairwiseDiffusionForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    D_(this->coeffs().getScalar("D")),
    h_(this->coeffs().getScalar("h"))
{}

template<class CloudType>
PairwiseDiffusionForce<CloudType>::PairwiseDiffusionForce
(
    const PairwiseDiffusionForce& dff
)
:
    ParticleForce<CloudType>(dff),
    D_(dff.D_),
    h_(dff.h_)
{}
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class CloudType>
Foam::forceSuSp Foam::PairwiseDiffusionForce<CloudType>::calcNonCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    const int myRank = Pstream::myProcNo();
    const int nProcs = Pstream::nProcs();

    CloudType& cloud = const_cast<CloudType&>(this->owner());
    const auto& mesh = cloud.mesh();
    auto& cellOccupancy = cloud.cellOccupancy();

    const label cellI = p.cell();
    const point& pCellCentre = mesh.C()[cellI];
    const scalar sigma = h_ / 2.0;

    vector totalForce = vector::zero;
    scalar weightSum  = 0.0;

    // ---- static cache updated once per timestep ----
    static label lastTimeIndex = -1;
    static List<point> remoteCenters;
    static List<label> remoteCounts;

    const Time& runTime = mesh.time();
    const label timeIndex = runTime.timeIndex();

    if (timeIndex != lastTimeIndex)
    {
        lastTimeIndex = timeIndex;

        remoteCenters.clear();
        remoteCounts.clear();

        // Build local data
        List<point> localCenters;
        List<label> localCounts;

        forAll(mesh.C(), cI)
        {
            if (cellOccupancy[cI].size() > 0)
            {
                localCenters.append(mesh.C()[cI]);
                localCounts.append(cellOccupancy[cI].size());
            }
        }

        // Use non-blocking PstreamBuffers + UOPstream/UIPstream to avoid
        // pairwise deadlock when running on many processors. This mirrors
        // the pattern used elsewhere (e.g. RecycleInteraction).
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send local lists to all other procs using UOPstream wrappers
        PtrList<UOPstream> UOPstreamPtrs(Pstream::nProcs());

        for (int procI = 0; procI < nProcs; ++procI)
        {
            if (procI == myRank) continue;

            auto* osptr = UOPstreamPtrs.get(procI);
            if (!osptr)
            {
                osptr = new UOPstream(procI, pBufs);
                UOPstreamPtrs.set(procI, osptr);
            }

            (*osptr) << localCenters << localCounts;

            Pout << "Proc " << myRank
                << " sent " << localCenters.size()
                << " cells to proc " << procI
                << " at time " << runTime.timeName() << endl;
        }

        // Ensure all sends are finished and incoming buffers are available
        pBufs.finishedSends();

        // Retrieve data from any procs that sent to us
        for (const int proci : pBufs.allProcs())
        {
            if (pBufs.recvDataCount(proci))
            {
                UIPstream is(proci, pBufs);
                List<point> recCenters;
                List<label> recCounts;
                is >> recCenters >> recCounts;

                Pout << "Proc " << myRank
                    << " received " << recCenters.size()
                    << " cells from proc " << proci
                    << " at time " << runTime.timeName() << endl;

                forAll(recCenters, i)
                {
                    remoteCenters.append(recCenters[i]);
                    remoteCounts.append(recCounts[i]);
                }
            }
        }
    }

    // ---- remote contributions ----
    forAll(remoteCenters, j)
    {
        const point& otherCentre = remoteCenters[j];
        const label nParcels     = remoteCounts[j];

        if (magSqr(pCellCentre - otherCentre) > sqr(h_)) continue;

        const vector dr = p.position() - otherCentre;
        const scalar distSqr = magSqr(dr);
        if (distSqr > sqr(h_)) continue;
        if (mag(dr) < VSMALL) continue;

        const scalar w =
            1.0 / sqrt(2.0 * constant::mathematical::pi * sqr(sigma))
          * exp(-0.5 * sqr(mag(dr)/sigma));

        totalForce += nParcels * w * (dr / mag(dr));
        weightSum  += nParcels * w;
    }

    // ---- local contributions (exact parcel interactions) ----
    forAll(mesh.C(), otherCellI)
    {
        const point& otherCentre = mesh.C()[otherCellI];
        if (magSqr(pCellCentre - otherCentre) > sqr(h_)) continue;

        const auto& parcelList = cellOccupancy[otherCellI];

        for (const auto* qPtr : parcelList)
        {
            if (qPtr == &p) continue;

            const vector dr = p.position() - qPtr->position();
            const scalar distSqr = magSqr(dr);
            if (distSqr > sqr(h_)) continue;
            if (mag(dr) < VSMALL) continue;

            const scalar w =
                1.0 / sqrt(2.0 * constant::mathematical::pi * sqr(sigma))
              * exp(-0.5 * sqr(mag(dr)/sigma));

            totalForce += w * (dr / mag(dr));
            weightSum  += w;
        }
    }

    // ---- final averaged force ----
    vector averagedForce =
        (weightSum > VSMALL ? D_ * totalForce / weightSum : vector::zero);

    return Foam::forceSuSp(averagedForce, 0.0);
}



// ************************************************************************* //
