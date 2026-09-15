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

#include "TwoPhaseThermophoreticForce.H"
#include "fvcGrad.H"
#include "mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TwoPhaseThermophoreticForce<CloudType>::TwoPhaseThermophoreticForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& forceType
)
:
    ParticleForce<CloudType>(owner, mesh, dict, forceType, true),
    T1Name_(this->coeffs().template getOrDefault<word>("T1", "T.water")),
    T2Name_(this->coeffs().template getOrDefault<word>("T2", "T.air")),
    Kp_(readScalar(this->coeffs().lookup("Kp"))),
    Kc1_(readScalar(this->coeffs().lookup("Kc1"))),
    Kc2_(readScalar(this->coeffs().lookup("Kc2"))),
    lambda1_(readScalar(this->coeffs().lookup("lambda1"))),
    lambda2_(readScalar(this->coeffs().lookup("lambda2"))),
    gradT1InterpPtr_(nullptr),
    gradT2InterpPtr_(nullptr)
{}


template<class CloudType>
Foam::TwoPhaseThermophoreticForce<CloudType>::TwoPhaseThermophoreticForce
(
    const TwoPhaseThermophoreticForce& pgf
)
:
    ParticleForce<CloudType>(pgf),
    T1Name_(pgf.T1Name_),
    T2Name_(pgf.T2Name_),
    Kp_(pgf.Kp_),
    Kc1_(pgf.Kc1_),
    Kc2_(pgf.Kc2_),
    lambda1_(pgf.lambda1_),
    lambda2_(pgf.lambda2_),
    gradT1InterpPtr_(nullptr),
    gradT2InterpPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::TwoPhaseThermophoreticForce<CloudType>::cacheFields
(
    const bool store
)
{
    if (store)
    {
        {
            static word fName("gradT1");

            if (!this->mesh().template foundObject<volVectorField>(fName))
            {
                const volScalarField& Tc = this->mesh().template
                    lookupObject<volScalarField>(T1Name_);

                volVectorField* gradTPtr =
                    new volVectorField(fName, fvc::grad(Tc));

                gradTPtr->store();
            }

            gradT1InterpPtr_.reset
            (
                interpolation<vector>::New
                (
                    this->owner().solution().interpolationSchemes(),
                    this->mesh().template lookupObject<volVectorField>(fName)
                ).ptr()
            );
        }

        {
            static word fName("gradT2");

            if (!this->mesh().template foundObject<volVectorField>(fName))
            {
                const volScalarField& Tc = this->mesh().template
                    lookupObject<volScalarField>(T2Name_);

                volVectorField* gradTPtr =
                    new volVectorField(fName, fvc::grad(Tc));

                gradTPtr->store();
            }

            gradT2InterpPtr_.reset
            (
                interpolation<vector>::New
                (
                    this->owner().solution().interpolationSchemes(),
                    this->mesh().template lookupObject<volVectorField>(fName)
                ).ptr()
            );
        }
    }
    else
    {
        gradT1InterpPtr_.clear();
        gradT2InterpPtr_.clear();

        for (const word& fName : {word("gradT1"), word("gradT2")})
        {
            if (this->mesh().template foundObject<volVectorField>(fName))
            {
                const_cast<volVectorField&>
                (
                    this->mesh().template lookupObject<volVectorField>(fName)
                ).checkOut();
            }
        }
    }
}


template<class CloudType>
Foam::forceSuSp Foam::TwoPhaseThermophoreticForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    const scalar Cs = 1.17;
    const scalar Ct = 2.18;
    const scalar Cm = 1.14;

    // Phase 1 (primary)
    const scalar Kn1 = 2.0*lambda1_/p.d();
    const scalar K1 = Kc1_/Kp_;
    const scalar Dt1 =
    (
        (-6.0*mathematical::pi*sqr(muc)*p.d()*Cs*(K1 + Ct*Kn1))
       /(td.rhoc()*(1.0 + 3.0*Cm*Kn1)*(1.0 + 2.0*K1 + 2.0*Ct*Kn1))
    );
    const vector gradT1 =
        gradT1Interp().interpolate(p.coordinates(), p.currentTetIndices());

    forceSuSp F1(Zero);
    F1.Su() = -1.0*Dt1*gradT1/td.Tc();

    // Phase 2 (secondary) - same correlation, secondary phase's own
    // viscosity/density/temperature-gradient/molecular mean free path/
    // thermal conductivity
    const scalar Kn2 = 2.0*lambda2_/p.d();
    const scalar K2 = Kc2_/Kp_;
    const scalar Dt2 =
    (
        (-6.0*mathematical::pi*sqr(td.mu2c())*p.d()*Cs*(K2 + Ct*Kn2))
       /(td.rho2c()*(1.0 + 3.0*Cm*Kn2)*(1.0 + 2.0*K2 + 2.0*Ct*Kn2))
    );
    const vector gradT2 =
        gradT2Interp().interpolate(p.coordinates(), p.currentTetIndices());

    forceSuSp F2(Zero);
    F2.Su() = -1.0*Dt2*gradT2/td.T2c();

    const scalar chi = p.chi();

    return (1.0 - chi)*F1 + chi*F2;
}


// ************************************************************************* //
