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

#include "FunctionYPlusDispersionRAS.H"
#include "constants.H"
#include "nearWallDist.H"
#include "wallFvPatch.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FunctionYPlusDispersionRAS<CloudType>::FunctionYPlusDispersionRAS
(
    const dictionary& dict,
    CloudType& owner
)
:
    DispersionRASModel<CloudType>(dict, owner),
    yP_
    (
        IOobject
        (
        "yP",
        this->owner_.mesh().time().timeName(),
        this->owner_.mesh(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
        this->owner_.mesh(),
        dimensionedScalar("yP",dimless,0.0)
    )
{}


template<class CloudType>
Foam::FunctionYPlusDispersionRAS<CloudType>::FunctionYPlusDispersionRAS
(
    const FunctionYPlusDispersionRAS<CloudType>& dm
)
:
    DispersionRASModel<CloudType>(dm),
    yP_(dm.yP_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FunctionYPlusDispersionRAS<CloudType>::~FunctionYPlusDispersionRAS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::vector Foam::FunctionYPlusDispersionRAS<CloudType>::update
(
    const scalar dt,
    const label celli,
    const vector& U,
    const scalar d,
    const scalar rho,
    const vector& Uc,
    const scalar rhoc,
    const scalar muc,
    const vector& vortc,
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

        // Calculate yPlus value

        const fvPatchList& patches = this->owner().mesh().boundary();
        const objectRegistry& obr = this->owner().mesh();
        const word turbName =
            IOobject::groupName
            (
                turbulenceModel::propertiesName,
                this->owner().U().group()
            );
        const turbulenceModel& turbModel =
            obr.lookupObject<turbulenceModel>(turbName);

        forAll(patches, patchi)
        {
            // Obtain the current patch
            const fvPatch& currPatch = patches[patchi];

            // Calculation
            if (isA<wallFvPatch>(currPatch))
            {
                const scalarField& y_ = turbModel.y()[patchi];
                const fvPatchVectorField& Uw =
                    turbModel.U().boundaryField()[patchi];
                const tmp<scalarField> tnuw = turbModel.nu(patchi);
                const scalarField& nuw = tnuw();

                // yPlus for the wallFvPatch
                fvPatchField<scalar>& yPpatch =
                    const_cast<fvPatchField<scalar>&>(yP_.boundaryField()[patchi]);
                yPpatch = (y_ * sqrt(nuw * mag(Uw.snGrad())) / nuw);

                if (yPpatch[patchi] < 100.0)
                {
                    scalar f_x;
                    scalar f_y;
                    scalar f_z;

                    scalar dir_x;
                    scalar dir_y;
                    scalar dir_z;

                    const scalar sigma = sqrt(2 * k / 3.0);
                    const scalar theta = rnd.sample01<scalar>() * twoPi;
                    const scalar u = 2 * rnd.sample01<scalar>() - 1;
                    const scalar a = sqrt(1 - sqr(u));

                    f_x = 1 + 0.285 * (yPpatch[patchi] + 6) * exp(-0.455 * pow(yPpatch[patchi] + 6, 0.53));
                    f_y = 1 - exp(-0.02 * yPpatch[patchi]);
                    f_z = sqrt(3 - sqr(f_x) - sqr(f_y));

                    // Update the turbulence fluctuating velocity for particle
                    dir_x = f_x * a * cos(theta);
                    dir_y = f_y * a * sin(theta);
                    dir_z = f_z * u;

                    const vector dir(dir_x, dir_y, dir_z);

                    UTurb = sigma * dir;
                }
                else
                {
                    const scalar sigma = sqrt(2 * k / 3.0);

                    // Calculate a random exit direction `dir` distributed uniformly
                    // in spherical coordinates
                    const scalar theta = rnd.sample01<scalar>() * twoPi;
                    const scalar u = 2 * rnd.sample01<scalar>() - 1;

                    const scalar a = sqrt(1 - sqr(u));
                    const vector dir(a * cos(theta), a * sin(theta), u);

                    UTurb = sigma * mag(rnd.GaussNormal<scalar>()) * dir;
                }
            }
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
