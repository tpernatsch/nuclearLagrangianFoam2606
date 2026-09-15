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

#include "TwoPhasePressureGradientForce.H"
#include "fvcDdt.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TwoPhasePressureGradientForce<CloudType>::TwoPhasePressureGradientForce
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
    DU1cDtInterpPtr_(nullptr),
    DU2cDtInterpPtr_(nullptr)
{}


template<class CloudType>
Foam::TwoPhasePressureGradientForce<CloudType>::TwoPhasePressureGradientForce
(
    const TwoPhasePressureGradientForce& pgf
)
:
    ParticleForce<CloudType>(pgf),
    U1Name_(pgf.U1Name_),
    U2Name_(pgf.U2Name_),
    DU1cDtInterpPtr_(nullptr),
    DU2cDtInterpPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::TwoPhasePressureGradientForce<CloudType>::cacheFields
(
    const bool store
)
{
    if (store)
    {
        {
            static word resultName("DU1cDt");

            volVectorField* resultPtr =
                this->mesh().template getObjectPtr<volVectorField>(resultName);

            if (!resultPtr)
            {
                const volVectorField& Uc = this->mesh().template
                    lookupObject<volVectorField>(U1Name_);

                resultPtr = new volVectorField
                (
                    resultName,
                    fvc::ddt(Uc) + (Uc & fvc::grad(Uc))
                );

                resultPtr->store();
            }

            DU1cDtInterpPtr_.reset
            (
                interpolation<vector>::New
                (
                    this->owner().solution().interpolationSchemes(),
                    *resultPtr
                ).ptr()
            );
        }

        {
            static word resultName("DU2cDt");

            volVectorField* resultPtr =
                this->mesh().template getObjectPtr<volVectorField>(resultName);

            if (!resultPtr)
            {
                const volVectorField& Uc = this->mesh().template
                    lookupObject<volVectorField>(U2Name_);

                resultPtr = new volVectorField
                (
                    resultName,
                    fvc::ddt(Uc) + (Uc & fvc::grad(Uc))
                );

                resultPtr->store();
            }

            DU2cDtInterpPtr_.reset
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
        DU1cDtInterpPtr_.clear();
        DU2cDtInterpPtr_.clear();

        for (const word& resultName : {word("DU1cDt"), word("DU2cDt")})
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
Foam::forceSuSp Foam::TwoPhasePressureGradientForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    const vector DU1cDt =
        DU1cDtInterp().interpolate(p.coordinates(), p.currentTetIndices());
    const vector DU2cDt =
        DU2cDtInterp().interpolate(p.coordinates(), p.currentTetIndices());

    forceSuSp F1(Zero);
    F1.Su() = mass*td.rhoc()/p.rho()*DU1cDt;

    forceSuSp F2(Zero);
    F2.Su() = mass*td.rho2c()/p.rho()*DU2cDt;

    const scalar chi = p.chi();

    return (1.0 - chi)*F1 + chi*F2;
}


template<class CloudType>
Foam::scalar Foam::TwoPhasePressureGradientForce<CloudType>::massAdd
(
    const typename CloudType::parcelType&,
    const typename CloudType::parcelType::trackingData&,
    const scalar
) const
{
    return 0.0;
}


// ************************************************************************* //
