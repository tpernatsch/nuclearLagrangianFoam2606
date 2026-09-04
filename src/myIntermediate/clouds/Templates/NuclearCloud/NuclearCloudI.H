/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
using namespace Foam::constant;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::NuclearCloud<CloudType>&
Foam::NuclearCloud<CloudType>::cloudCopy() const
{
    return *cloudCopyPtr_;
}


template<class CloudType>
inline const typename CloudType::particleType::constantProperties&
Foam::NuclearCloud<CloudType>::constProps() const
{
    return constProps_;
}


template<class CloudType>
inline typename CloudType::particleType::constantProperties&
Foam::NuclearCloud<CloudType>::constProps()
{
    return constProps_;
}


template<class CloudType>
inline const Foam::volVectorField& Foam::NuclearCloud<CloudType>::vort() const
{
    return vort_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::NuclearCloud<CloudType>::T() const
{
    return T_;
}

template<class CloudType>
inline const Foam::HeatTransferModel<Foam::NuclearCloud<CloudType>>&
Foam::NuclearCloud<CloudType>::heatTransfer() const
{
    return *heatTransferModel_;
}


template<class CloudType>
inline const Foam::DecayHeatModel<Foam::NuclearCloud<CloudType>>&
Foam::NuclearCloud<CloudType>::decayHeat() const
{
    return *decayHeatModel_;
}


template<class CloudType>
inline const Foam::integrationScheme&
Foam::NuclearCloud<CloudType>::TIntegrator() const
{
    return *TIntegrator_;
}


template<class CloudType>
inline bool Foam::NuclearCloud<CloudType>::radiation() const
{
    return radiation_;
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::NuclearCloud<CloudType>::radAreaP()
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return *radAreaP_;
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::NuclearCloud<CloudType>::radAreaP() const
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return *radAreaP_;
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::NuclearCloud<CloudType>::radT4()
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return *radT4_;
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::NuclearCloud<CloudType>::radT4() const
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return *radT4_;
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::NuclearCloud<CloudType>::radAreaPT4()
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return *radAreaPT4_;
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::NuclearCloud<CloudType>::radAreaPT4() const
{
    if (!radiation_)
    {
        FatalErrorInFunction
            << "Radiation field requested, but radiation model not active"
            << abort(FatalError);
    }

    return *radAreaPT4_;
}


template<class CloudType>
inline void Foam::NuclearCloud<CloudType>::transferToCarrier
(
    const parcelType& p,
    const typename parcelType::trackingData& td
)
{
    CloudType::transferToCarrier(p, td);

    hsTrans()[p.cell()] += p.nParticle()*p.mass()*p.hs();
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::NuclearCloud<CloudType>::hsTrans()
{
    return *hsTrans_;
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::NuclearCloud<CloudType>::hsTrans() const
{
    return *hsTrans_;
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::NuclearCloud<CloudType>::hsCoeff()
{
    return *hsCoeff_;
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::NuclearCloud<CloudType>::hsCoeff() const
{
    return *hsCoeff_;
}

template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::NuclearCloud<CloudType>::hsCoeffTemp()
{
    return *hsCoeffTemp_;
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::NuclearCloud<CloudType>::hsCoeffTemp() const
{
    return *hsCoeffTemp_;
}


template<class CloudType>
inline Foam::tmp<Foam::fvScalarMatrix>
Foam::NuclearCloud<CloudType>::Sh(volScalarField& T, const dimensionedScalar& rho,  const dimensionedScalar& Cp) const
{
    const dimensionedScalar rhoCp = rho*Cp;

    DebugInfo
        << "hsTrans min/max = " << min(hsTrans()).value() << ", "
        << max(hsTrans()).value() << nl
        << "hsCoeff min/max = " << min(hsCoeff()).value() << ", "
        << max(hsCoeff()).value() << endl;

    const volScalarField::Internal
                Vdt(this->mesh().V()*this->db().time().deltaT());

    if (this->solution().coupledT())
    {
        if (this->solution().semiImplicit("T"))
        {
            return
              //  hsTrans()/(Cp)/(this->db().time().deltaT())
              //- fvm::SuSp(hsCoeff()/(Vdt*rhoCp), T)
              //+ hsCoeff()/(Vdt*rhoCp)*T;
                hsTrans()/(this->db().time().deltaT()*rhoCp) // Corrected source term to Temperature/Time
              - fvm::SuSp(hsCoeff()/rhoCp, T)                 // Corrected coefficient for fvm::SuSp to 1/Time
              + hsCoeff()/rhoCp*T;
        }
        else
        {
            //tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(T, dimTemperature*dimVolume/dimTime));
            tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(T, dimTemperature/dimTime));
            fvScalarMatrix& fvm = tfvm.ref();

            //fvm.source() = -hsTrans()/(Cp)/(this->db().time().deltaT());
            fvm.source() = -hsTrans()/(this->db().time().deltaT()*rhoCp);

            return tfvm;
        }
    }

    return tmp<fvScalarMatrix>::New(T, dimTemperature/dimTime);
}


template<class CloudType>
inline Foam::tmp<Foam::volScalarField> Foam::NuclearCloud<CloudType>::Ep() const
{
    tmp<volScalarField> tEp
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":radiation:Ep",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            this->mesh(),
            dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero)
        )
    );

    if (radiation_)
    {
        scalarField& Ep = tEp.ref().primitiveFieldRef();
        const scalar dt = this->db().time().deltaTValue();
        const scalarField& V = this->mesh().V();
        const scalar epsilon = constProps_.epsilon0();
        const scalarField& sumAreaPT4 = radAreaPT4_->field();

        Ep = sumAreaPT4*epsilon*physicoChemical::sigma.value()/V/dt;
    }

    return tEp;
}


template<class CloudType>
inline Foam::tmp<Foam::volScalarField> Foam::NuclearCloud<CloudType>::ap() const
{
    tmp<volScalarField> tap
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":radiation:ap",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            this->mesh(),
            dimensionedScalar(dimless/dimLength, Zero)
        )
    );

    if (radiation_)
    {
        scalarField& ap = tap.ref().primitiveFieldRef();
        const scalar dt = this->db().time().deltaTValue();
        const scalarField& V = this->mesh().V();
        const scalar epsilon = constProps_.epsilon0();
        const scalarField& sumAreaP = radAreaP_->field();

        ap = sumAreaP*epsilon/V/dt;
    }

    return tap;
}


template<class CloudType>
inline Foam::tmp<Foam::volScalarField>
Foam::NuclearCloud<CloudType>::sigmap() const
{
    tmp<volScalarField> tsigmap
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":radiation:sigmap",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            this->mesh(),
            dimensionedScalar(dimless/dimLength, Zero)
        )
    );

    if (radiation_)
    {
        scalarField& sigmap = tsigmap.ref().primitiveFieldRef();
        const scalar dt = this->db().time().deltaTValue();
        const scalarField& V = this->mesh().V();
        const scalar epsilon = constProps_.epsilon0();
        const scalar f = constProps_.f0();
        const scalarField& sumAreaP = radAreaP_->field();

        sigmap = sumAreaP*(1.0 - f)*(1.0 - epsilon)/V/dt;
    }

    return tsigmap;
}


template<class CloudType>
inline Foam::scalar Foam::NuclearCloud<CloudType>::Tmax() const
{
    scalar val = -GREAT;
    bool nonEmpty = false;

    for (const parcelType& p : *this)
    {
        val = max(val, p.T());
        nonEmpty = true;
    }

    if (returnReduceOr(nonEmpty))
    {
        return returnReduce(val, maxOp<scalar>());
    }

    return 0;
}


template<class CloudType>
inline Foam::scalar Foam::NuclearCloud<CloudType>::Tmin() const
{
    scalar val = GREAT;
    bool nonEmpty = false;

    for (const parcelType& p : *this)
    {
        val = min(val, p.T());
        nonEmpty = true;
    }

    if (returnReduceOr(nonEmpty))
    {
        return returnReduce(val, minOp<scalar>());
    }

    return 0;
}


// ************************************************************************* //
