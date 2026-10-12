/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "AnisotropicStochasticDispersionRAS.H"
#include "constants.H"
#include "turbulenceModel.H"
#include "wallFvPatch.H"
#include "wallDistAddressing.H"
#include "mapDistribute.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{
    //- An arbitrary unit vector perpendicular to n, used as a fallback
    //  streamwise direction where the local mean flow has (almost) no
    //  tangential component, e.g. at a stagnation point
    inline vector anisotropicStochasticDispersionRAS_arbitraryPerp
    (
        const vector& n
    )
    {
        vector a(1, 0, 0);
        if (mag(n.x()) > mag(n.y()) && mag(n.x()) > mag(n.z()))
        {
            a = vector(0, 1, 0);
        }

        vector t(a - (a & n)*n);
        return t/(mag(t) + ROOTVSMALL);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::AnisotropicStochasticDispersionRAS<CloudType>::
AnisotropicStochasticDispersionRAS
(
    const dictionary& dict,
    CloudType& owner
)
:
    DispersionRASModel<CloudType>(dict, owner),
    yPlusThreshold_
    (
        this->coeffDict().template getOrDefault<scalar>("yPlusThreshold", 100)
    ),
    yPtr_(nullptr),
    uTauPtr_(nullptr),
    nuPtr_(nullptr),
    nWallPtr_(nullptr)
{}


template<class CloudType>
Foam::AnisotropicStochasticDispersionRAS<CloudType>::
AnisotropicStochasticDispersionRAS
(
    const AnisotropicStochasticDispersionRAS<CloudType>& dm
)
:
    DispersionRASModel<CloudType>(dm),
    yPlusThreshold_(dm.yPlusThreshold_),
    yPtr_(nullptr),
    uTauPtr_(nullptr),
    nuPtr_(nullptr),
    nWallPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::AnisotropicStochasticDispersionRAS<CloudType>::
~AnisotropicStochasticDispersionRAS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::AnisotropicStochasticDispersionRAS<CloudType>::cacheFields
(
    const bool store
)
{
    DispersionRASModel<CloudType>::cacheFields(store);

    if (!store)
    {
        yPtr_.clear();
        uTauPtr_.clear();
        nuPtr_.clear();
        nWallPtr_.clear();
        return;
    }

    const fvMesh& mesh = this->owner().mesh();

    const objectRegistry& obr = mesh;
    const word turbName =
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            this->owner().U().group()
        );

    const turbulenceModel* turb = obr.findObject<turbulenceModel>(turbName);

    if (!turb)
    {
        FatalErrorInFunction
            << "Turbulence model not found in mesh database" << nl
            << "Database objects include: " << obr.sortedToc()
            << abort(FatalError);
    }

    nuPtr_.reset(turb->nu().ptr());

    // Wall-patch friction velocity and wall-normal (into the fluid),
    // left as the default (zero) value on non-wall patches and in the
    // domain interior - filled in below via wallDistAddressing::map()
    auto uTauPtr = autoPtr<volScalarField>::New
    (
        IOobject
        (
            "AnisotropicStochasticDispersionRAS:uTau",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimVelocity, Zero)
    );

    auto nWallPtr = autoPtr<volVectorField>::New
    (
        IOobject
        (
            "AnisotropicStochasticDispersionRAS:nWall",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedVector(dimless, Zero)
    );

    const volVectorField& U = turb->U();
    volScalarField::Boundary& uTauBf = uTauPtr->boundaryFieldRef();
    volVectorField::Boundary& nWallBf = nWallPtr->boundaryFieldRef();
    const fvPatchList& patches = mesh.boundary();

    forAll(patches, patchi)
    {
        if (isA<wallFvPatch>(patches[patchi]))
        {
            const fvPatchVectorField& Uw = U.boundaryField()[patchi];

            uTauBf[patchi] =
                sqrt(turb->nuEff(patchi)*mag(Uw.snGrad()) + VSMALL);

            nWallBf[patchi] = -patches[patchi].nf();
        }
    }

    // Transport the wall-patch values to the nearest cell/patch throughout
    // the domain
    const auto& wDist = wallDistAddressing::New(mesh);

    yPtr_.reset(new volScalarField(wDist.y()));

    wDist.map(uTauPtr(), mapDistribute::transform());
    wDist.map(nWallPtr(), mapDistribute::transform());

    uTauPtr_ = std::move(uTauPtr);
    nWallPtr_ = std::move(nWallPtr);
}


template<class CloudType>
Foam::vector Foam::AnisotropicStochasticDispersionRAS<CloudType>::update
(
    const scalar dt,
    const label celli,
    const vector& U,
    const scalar d,
    const scalar rho,
    const vector& Uc,
    const scalar rhoc,
    const scalar muc,
    vector& UTurb,
    scalar& tTurb,
    scalar& tTurbLoc
)
{
    Random& rnd = this->owner().rndGen();

    const scalar cps = 0.16432;

    const scalar k = this->kPtr_->primitiveField()[celli];
    const scalar epsilon =
        this->epsilonPtr_->primitiveField()[celli] + ROOTVSMALL;

    const scalar UrelMag = mag(U - Uc - UTurb);

    tTurbLoc =
        min(k/epsilon, cps*pow(k, 1.5)/epsilon/(UrelMag + SMALL));

    // Parcel is perturbed by the turbulence
    if (dt < tTurbLoc)
    {
        tTurb += dt;

        if (tTurb > tTurbLoc)
        {
            tTurb = 0;

            const scalar y = yPtr_->primitiveField()[celli];
            const scalar uTau = uTauPtr_->primitiveField()[celli];
            const scalar nu = nuPtr_->primitiveField()[celli] + ROOTVSMALL;

            const scalar yPlus = y*uTau/nu;

            if (yPlus < yPlusThreshold_ && y < 0.5*GREAT)
            {
                // Near-wall: anisotropic draw in a local wall-oriented
                // frame (Dreeben & Pope, 1997)

                const vector& n = nWallPtr_->primitiveField()[celli];

                vector streamwise(Uc - (Uc & n)*n);
                const scalar magStreamwise = mag(streamwise);

                if (magStreamwise > 1e-3*max(mag(Uc), SMALL))
                {
                    streamwise /= magStreamwise;
                }
                else
                {
                    // Degenerate case (e.g. stagnation point): mean flow
                    // has (almost) no tangential component, so pick an
                    // arbitrary stable tangential direction instead
                    streamwise =
                        anisotropicStochasticDispersionRAS_arbitraryPerp(n);
                }

                const vector spanwise(n ^ streamwise);

                const scalar uPlusRMS =
                    0.4*yPlus/(1 + 0.0239*pow(yPlus, 1.496));
                const scalar vPlusRMS =
                    0.116*sqr(yPlus)
                   /(1 + 0.0203*yPlus + 0.0014*pow(yPlus, 2.421));
                const scalar wPlusRMS =
                    0.19*yPlus/(1 + 0.0361*pow(yPlus, 1.322));

                UTurb =
                    uPlusRMS*uTau*rnd.GaussNormal<scalar>()*streamwise
                  + vPlusRMS*uTau*rnd.GaussNormal<scalar>()*n
                  + wPlusRMS*uTau*rnd.GaussNormal<scalar>()*spanwise;
            }
            else
            {
                // Bulk: isotropic, as StochasticDispersionRAS

                const scalar sigma = sqrt(2*k/3.0);

                // Calculate a random direction dir distributed uniformly
                // in spherical coordinates

                const scalar theta = rnd.sample01<scalar>()*twoPi;
                const scalar uDir = 2*rnd.sample01<scalar>() - 1;

                const scalar a = sqrt(1 - sqr(uDir));
                const vector dir(a*cos(theta), a*sin(theta), uDir);

                UTurb = sigma*mag(rnd.GaussNormal<scalar>())*dir;
            }
        }
    }
    else
    {
        tTurb = GREAT;
        UTurb = Zero;
    }

    return Uc + UTurb;
}


// ************************************************************************* //
