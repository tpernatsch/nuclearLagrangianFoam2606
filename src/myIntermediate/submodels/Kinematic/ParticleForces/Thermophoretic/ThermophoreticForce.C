/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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
#include "ThermophoreticForce.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //



template<class CloudType>
Foam::ThermophoreticForce<CloudType>::ThermophoreticForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& forceType
)
:
    ParticleForce<CloudType>(owner, mesh, dict, forceType, true),
    TName_(this->coeffs().template lookupOrDefault<word>("T", "T")),
    Kp_(readScalar(this->coeffs().lookup("Kp"))),
    Kc_(readScalar(this->coeffs().lookup("Kc"))),
    lambda_(readScalar(this->coeffs().lookup("lambda"))),
    gradTInterpPtr_(nullptr)
{}


template<class CloudType>
Foam::ThermophoreticForce<CloudType>::ThermophoreticForce
(
    const ThermophoreticForce& pgf
)
:
    ParticleForce<CloudType>(pgf),
    TName_(pgf.TName_),
    Kp_(pgf.Kp_),
    Kc_(pgf.Kc_),
    lambda_(pgf.lambda_),
    gradTInterpPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermophoreticForce<CloudType>::~ThermophoreticForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ThermophoreticForce<CloudType>::cacheFields(const bool store)
{
    static word fName("gradT");

    bool fieldExists = this->mesh().template foundObject<volVectorField>(fName);

    if (store)
    {
        if (!fieldExists)
        {
            const volScalarField& Tc = this->mesh().template
                lookupObject<volScalarField>(TName_);

            volVectorField* gradTPtr = new volVectorField
            (
                fName,
                (fvc::grad(Tc))
            );

            gradTPtr->store();
        }

        const volVectorField& gradT = this->mesh().template
            lookupObject<volVectorField>(fName);

        gradTInterpPtr_.reset
        (
            interpolation<vector>::New
            (
                this->owner().solution().interpolationSchemes(),
                gradT
            ).ptr()
        );
    }
    else
    {
        gradTInterpPtr_.clear();

        if (fieldExists)
        {
            const volVectorField& gradT = this->mesh().template
                lookupObject<volVectorField>(fName);

            const_cast<volVectorField&>(gradT).checkOut();
        }
    }
   
}


template<class CloudType>
Foam::forceSuSp Foam::ThermophoreticForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(Zero, 0.0);
    const scalar Cs=1.17;
    const scalar Ct=2.18;
    const scalar Cm=1.14;
    const scalar Kn=2*lambda_/p.d();
    const scalar K=Kc_/Kp_;
    
    const scalar Dt=((-6.0*mathematical::pi*sqr(muc)*p.d()*Cs*(K+Ct*Kn))/(td.rhoc()*(1.0+3.0*Cm*Kn)*(1.0+2.0*K+2.0*Ct*Kn)));
    vector gradT =
        gradTInterp().interpolate(p.coordinates(), p.currentTetIndices());

    value.Su() = -1.0*Dt*gradT/td.Tc();
  
    return value;
     
}

// ************************************************************************* //
