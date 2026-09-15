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

\*---------------------------------------------------------------------------*/

#include "TwoPhaseGradientDispersionRAS.H"
#include "demandDrivenData.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TwoPhaseGradientDispersionRAS<CloudType>::TwoPhaseGradientDispersionRAS
(
    const dictionary& dict,
    CloudType& owner
)
:
    TwoPhaseDispersionRASModel<CloudType>(dict, owner),
    gradk2Ptr_(nullptr),
    ownGradK2_(false)
{}


template<class CloudType>
Foam::TwoPhaseGradientDispersionRAS<CloudType>::TwoPhaseGradientDispersionRAS
(
    const TwoPhaseGradientDispersionRAS<CloudType>& dm
)
:
    TwoPhaseDispersionRASModel<CloudType>(dm),
    gradk2Ptr_(dm.gradk2Ptr_),
    ownGradK2_(dm.ownGradK2_)
{
    dm.ownGradK2_ = false;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TwoPhaseGradientDispersionRAS<CloudType>::~TwoPhaseGradientDispersionRAS()
{
    cacheFields(false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::TwoPhaseGradientDispersionRAS<CloudType>::cacheFields
(
    const bool store
)
{
    if (store)
    {
        gradk2Ptr_ = fvc::grad(this->k2Field()).ptr();
        ownGradK2_ = true;
    }
    else
    {
        if (ownGradK2_)
        {
            deleteDemandDrivenData(gradk2Ptr_);
            gradk2Ptr_ = nullptr;
            ownGradK2_ = false;
        }
    }
}


template<class CloudType>
Foam::vector Foam::TwoPhaseGradientDispersionRAS<CloudType>::update
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

    const scalar k = this->k2Field().primitiveField()[celli];
    const scalar epsilon =
        this->epsilon2Field().primitiveField()[celli] + ROOTVSMALL;
    const vector& gradk = this->gradk2Ptr_->primitiveField()[celli];

    const scalar UrelMag = mag(U - Uc - UTurb);

    tTurbLoc =
        min(k/epsilon, cps*pow(k, 1.5)/epsilon/(UrelMag + SMALL));


    // Parcel is perturbed by the turbulence
    if (dt < tTurbLoc)
    {
        tTurb += dt;

        if (tTurb > tTurbLoc)
        {
            tTurb = 0.0;

            const scalar sigma = sqrt(2.0*k/3.0);
            const vector dir = -gradk/(mag(gradk) + SMALL);

            scalar fac = 0.0;

            // In 2D calculations the -grad(k2) is always
            // away from the axis of symmetry
            // This creates a 'hole' in the spray and to
            // prevent this we let fac be both negative/positive
            if (this->owner().mesh().nSolutionD() == 2)
            {
                fac = rnd.GaussNormal<scalar>();
            }
            else
            {
                fac = mag(rnd.GaussNormal<scalar>());
            }

            UTurb = sigma*fac*dir;
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
