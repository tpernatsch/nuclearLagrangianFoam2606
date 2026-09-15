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

#include "TwoPhaseLiftForce.H"
#include "fvcCurl.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::TwoPhaseLiftForce<CloudType>::Cl
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const vector& curlUc,
    const scalar rhoc,
    const scalar Re,
    const scalar muc
) const
{
    // dummy
    return 0.0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TwoPhaseLiftForce<CloudType>::TwoPhaseLiftForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& forceType
)
:
    ParticleForce<CloudType>(owner, mesh, dict, forceType, true),
    U1Name_(this->coeffs().template getOrDefault<word>("U1", "U.water")),
    U2Name_(this->coeffs().template getOrDefault<word>("U2", "U.air")),
    curlU1cInterpPtr_(nullptr),
    curlU2cInterpPtr_(nullptr)
{}


template<class CloudType>
Foam::TwoPhaseLiftForce<CloudType>::TwoPhaseLiftForce(const TwoPhaseLiftForce& lf)
:
    ParticleForce<CloudType>(lf),
    U1Name_(lf.U1Name_),
    U2Name_(lf.U2Name_),
    curlU1cInterpPtr_(nullptr),
    curlU2cInterpPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::TwoPhaseLiftForce<CloudType>::cacheFields(const bool store)
{
    if (store)
    {
        {
            static word resultName("curlU1c");

            volVectorField* resultPtr =
                this->mesh().template getObjectPtr<volVectorField>(resultName);

            if (!resultPtr)
            {
                const volVectorField& Uc = this->mesh().template
                    lookupObject<volVectorField>(U1Name_);

                resultPtr = new volVectorField(resultName, fvc::curl(Uc));
                resultPtr->store();
            }

            curlU1cInterpPtr_.reset
            (
                interpolation<vector>::New
                (
                    this->owner().solution().interpolationSchemes(),
                    *resultPtr
                ).ptr()
            );
        }

        {
            static word resultName("curlU2c");

            volVectorField* resultPtr =
                this->mesh().template getObjectPtr<volVectorField>(resultName);

            if (!resultPtr)
            {
                const volVectorField& Uc = this->mesh().template
                    lookupObject<volVectorField>(U2Name_);

                resultPtr = new volVectorField(resultName, fvc::curl(Uc));
                resultPtr->store();
            }

            curlU2cInterpPtr_.reset
            (
                interpolation<vector>::New
                (
                    this->owner().solution().interpolationSchemes(),
                    *resultPtr
                ).ptr()
            );
        }
    }
    else
    {
        curlU1cInterpPtr_.clear();
        curlU2cInterpPtr_.clear();

        for (const word& resultName : {word("curlU1c"), word("curlU2c")})
        {
            volVectorField* resultPtr =
                this->mesh().template getObjectPtr<volVectorField>(resultName);

            if (resultPtr)
            {
                resultPtr->checkOut();
            }
        }
    }
}


template<class CloudType>
Foam::forceSuSp Foam::TwoPhaseLiftForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    const vector curlU1c =
        curlU1cInterp().interpolate(p.coordinates(), p.currentTetIndices());
    const vector curlU2c =
        curlU2cInterp().interpolate(p.coordinates(), p.currentTetIndices());

    // Phase 2 slip Reynolds number, kept on the same rho0/mu0 footing as
    // Re1 (the caller-supplied Re) - see the accompanying report
    const scalar rho0 = this->owner().constProps().rho0();
    const scalar Re2 = max(p.Re(rho0, p.U(), td.U2c(), p.d(), muc), SMALL);

    const scalar Cl1 = this->Cl(p, td, curlU1c, td.rhoc(), Re, muc);
    const scalar Cl2 = this->Cl(p, td, curlU2c, td.rho2c(), Re2, muc);

    forceSuSp F1(Zero);
    F1.Su() = mass/p.rho()*td.rhoc()*Cl1*((td.Uc() - p.U())^curlU1c);

    forceSuSp F2(Zero);
    F2.Su() = mass/p.rho()*td.rho2c()*Cl2*((td.U2c() - p.U())^curlU2c);

    const scalar chi = p.chi();

    return (1.0 - chi)*F1 + chi*F2;
}


// ************************************************************************* //
