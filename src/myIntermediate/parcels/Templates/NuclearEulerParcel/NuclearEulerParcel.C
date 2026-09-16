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

#include "NuclearEulerParcel.H"
#include "physicoChemicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::NuclearEulerParcel<ParcelType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    ParcelType::setCellValues(cloud, td);

    const tetIndices tetIs = this->currentTetIndices();

    // Void fractions - clamped into [0,1] since the underlying cell/point
    // interpolation scheme is not guaranteed to be strictly bounded
    td.alpha1c() =
        clamp(td.alpha1Interp().interpolate(this->coordinates(), tetIs), 0, 1);
    td.alpha2c() =
        clamp(td.alpha2Interp().interpolate(this->coordinates(), tetIs), 0, 1);

    td.U2c() = td.U2Interp().interpolate(this->coordinates(), tetIs);
    td.rho2c() = td.rho2Interp().interpolate(this->coordinates(), tetIs);
    td.mu2c() = td.mu2Interp().interpolate(this->coordinates(), tetIs);

    td.T2c() = td.T2Interp().interpolate(this->coordinates(), tetIs);
    if (td.T2c() < cloud.constProps().TMin())
    {
        if (debug)
        {
            WarningInFunction
                << "Limiting observed secondary-phase temperature in cell "
                << this->cell() << " to " << cloud.constProps().TMin()
                << nl << endl;
        }

        td.T2c() = cloud.constProps().TMin();
    }

    // Turbulence closure fields - floored at zero (a negative
    // interpolated k/epsilon has no physical meaning)
    td.k1c() =
        max(td.k1Interp().interpolate(this->coordinates(), tetIs), 0);
    td.epsilon1c() =
        max(td.epsilon1Interp().interpolate(this->coordinates(), tetIs), 0);
    td.k2c() =
        max(td.k2Interp().interpolate(this->coordinates(), tetIs), 0);
    td.epsilon2c() =
        max(td.epsilon2Interp().interpolate(this->coordinates(), tetIs), 0);

    td.d2c() =
        max(td.d2Interp().interpolate(this->coordinates(), tetIs), 0);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::NuclearEulerParcel<ParcelType>::updateChi
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    // One-way switch: once true, never reverts
    if (this->chi_)
    {
        return;
    }

    const scalar gamma = cloud.chiSwitchModel().gamma(*this, td);

    // Evaluated against the Eulerian step, not the Lagrangian sub-step
    const scalar deltaT = cloud.mesh().time().deltaTValue();
    const scalar pSwitch = 1.0 - exp(-gamma*deltaT);

    if (cloud.rndGen().template sample01<scalar>() < pSwitch)
    {
        this->chi_ = true;
    }
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::NuclearEulerParcel<ParcelType>::calcDispersion
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    // Primary phase's own dispersion is always applied, unchanged
    ParcelType::calcDispersion(cloud, td, dt);

    // Secondary phase's dispersion, only when switched and RAS
    if (this->chi_ && cloud.k2IsRAS())
    {
        td.U2c() = cloud.dispersion2().update
        (
            dt,
            this->cell(),
            this->U_,
            this->d_,
            this->rho_,
            td.U2c(),
            td.rho2c(),
            td.mu2c(),
            this->UTurb_,
            this->tTurb_,
            this->tTurbLoc_
        );
    }
}


template<class ParcelType>
template<class TrackCloudType>
Foam::scalar Foam::NuclearEulerParcel<ParcelType>::calcHeatTransferTwoPhase
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const scalar Re1,
    const scalar Re2,
    const scalar Pr,
    const scalar kappa,
    const scalar NCpW,
    const scalar Sh,
    scalar& dhsTrans,
    scalar& Sph
)
{
    if (!cloud.heatTransfer().active())
    {
        return this->T_;
    }

    const scalar d = this->d();
    const scalar rho = cloud.constProps().rho0();
    const scalar As = this->areaS(d);
    const scalar V = this->volume(d);
    const scalar m = rho*V;
    const scalar Cp = cloud.constProps().Cp0();

    // Normalised void-fraction weights (alpha1 + alpha2 should already
    // sum to ~1, but the floor guards against interpolation error taking
    // both towards zero at the same point, e.g. right at a wall/free
    // surface where both fields are being blended)
    const scalar alphaSum =
        max(td.alpha1c() + td.alpha2c(), cloud.constProps().alphaSmall());
    const scalar w1 = td.alpha1c()/alphaSum;
    const scalar w2 = td.alpha2c()/alphaSum;

    // Per-phase convective heat transfer coefficients: the same Nu(Re,Pr)
    // correlation, evaluated once against each phase's own slip Re
    const scalar htc1 = cloud.heatTransfer().htc(d, Re1, Pr, kappa, NCpW);
    const scalar htc2 = cloud.heatTransfer().htc(d, Re2, Pr, kappa, NCpW);

    const scalar decHeat = cloud.decayHeat().decayPower(d);

    // Void-fraction-weighted effective heat transfer coefficient and an
    // htc-weighted mean carrier temperature. This is a two-resistance-
    // in-parallel closure that reduces exactly to plain
    // NuclearParcel::calcHeatTransfer() in the limit w2 -> 0.
    const scalar htcEff = w1*htc1 + w2*htc2;
    const scalar TcEff =
    (
        htcEff > ROOTVSMALL
      ? (w1*htc1*td.Tc() + w2*htc2*td.T2c())/htcEff
      : td.Tc()
    );

    // Calculate the integration coefficients
    const scalar bcp = htcEff*As/(m*Cp);
    const scalar acp = bcp*TcEff;

    scalar ancp = Sh;
    if (cloud.radiation())
    {
        const tetIndices tetIs = this->currentTetIndices();
        const scalar Gc = td.GInterp().interpolate(this->coordinates(), tetIs);
        const scalar sigma = physicoChemical::sigma.value();
        const scalar epsilon = cloud.constProps().epsilon0();

        ancp += As*epsilon*(Gc/4.0 - sigma*pow4(this->T_));
    }
    ancp += decHeat;
    ancp /= m*Cp;

    // Integrate to find the new parcel temperature
    const scalar deltaT =
        cloud.TIntegrator().delta(this->T_, dt, acp + ancp, bcp);
    const scalar deltaTncp = ancp*dt;
    const scalar deltaTcp = deltaT - deltaTncp;

    scalar Tnew = this->T_ + deltaT;
    Tnew = clamp(Tnew, cloud.constProps().TMin(), cloud.constProps().TMax());

    dhsTrans -= m*Cp*deltaTcp;
    Sph = dt*m*Cp*bcp;

    return Tnew;
}


// * * * * * * * * * * * * * * *  Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::NuclearEulerParcel<ParcelType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    this->updateChi(cloud, td);

    // Define local properties at beginning of time step
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar np0 = this->nParticle_;
    const scalar mass0 = this->mass();

    // Store T for consistent radiation source
    const scalar T0 = this->T_;


    // Calc surface values (based on the particle's own constant
    // properties, exactly as for plain NuclearParcel - carrier-phase
    // properties enter only through Re1/Re2/td.Tc()/td.T2c() below)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalar Ts, rhos, mus, Pr, kappas;
    this->calcSurfaceValues(cloud, td, this->T_, Ts, rhos, mus, Pr, kappas);

    // Slip Reynolds numbers against each Euler phase. Floored at SMALL,
    // as in plain NuclearParcel, to avoid a Stokes-equilibrated particle
    // underflowing Re to exactly 0 and faulting downstream in a
    // heat-transfer correlation.
    const scalar Re1 =
        max(this->Re(rhos, this->U_, td.Uc(), this->d_, mus), SMALL);
    const scalar Re2 =
        max(this->Re(rhos, this->U_, td.U2c(), this->d_, mus), SMALL);


    // Sources
    // ~~~~~~~

    // Explicit momentum source for particle
    vector Su = Zero;

    // Linearised momentum source coefficient
    scalar Spu = 0.0;

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = Zero;

    // Explicit enthalpy source for particle
    scalar Sh = 0.0;

    // Linearised source coefficient
    scalar Sph = 0.0;

    // Sensible enthalpy transfer from the particle to the carrier phase
    scalar dhsTrans = 0.0;


    // Heat transfer
    // ~~~~~~~~~~~~~

    // Sum Ni*Cpi*Wi of emission species
    scalar NCpW = 0.0;

    // Calculate new particle temperature from the void-fraction-weighted
    // two-phase closure
    this->T_ =
        this->calcHeatTransferTwoPhase
        (
            cloud,
            td,
            dt,
            Re1,
            Re2,
            Pr,
            kappas,
            NCpW,
            Sh,
            dhsTrans,
            Sph
        );


    // Motion
    // ~~~~~~

    // Calculate new particle velocity. Drag/buoyancy/other force models
    // act against the primary (tracking) phase only, exactly as for
    // plain NuclearParcel/KinematicParcel - unchanged.
    this->U_ =
        this->calcVelocity(cloud, td, dt, Re1, mus, mass0, Su, dUTrans, Spu);


    //  Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // One-way coupling only: these feed back onto the primary phase's
    // fields exactly like a plain NuclearParcel would, ready for a future
    // two-way-coupled fvOption/fvModel to consume. Nothing is currently
    // fed back onto the secondary phase.
    if (cloud.solution().coupledU() || cloud.solution().coupledT())
    {
        // Update momentum transfer
        cloud.UTrans()[this->cell()] += np0*dUTrans;

        // Update momentum transfer coefficient
        cloud.UCoeff()[this->cell()] += np0*Spu;

        // Update sensible enthalpy transfer
        cloud.hsTrans()[this->cell()] += np0*dhsTrans;

        // Update heat source coefficient
        cloud.hsCoeff()[this->cell()] += np0*Sph;

        // Update heat source coefficient multiplied by parcel temperature
        cloud.hsCoeffTemp()[this->cell()] += np0*Sph*this->T_;

        // Update radiation fields
        if (cloud.radiation())
        {
            const scalar ap = this->areaP();
            const scalar T4 = pow4(T0);
            cloud.radAreaP()[this->cell()] += dt*np0*ap;
            cloud.radT4()[this->cell()] += dt*np0*T4;
            cloud.radAreaPT4()[this->cell()] += dt*np0*ap*T4;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::NuclearEulerParcel<ParcelType>::NuclearEulerParcel
(
    const NuclearEulerParcel<ParcelType>& p
)
:
    ParcelType(p)
{}


template<class ParcelType>
Foam::NuclearEulerParcel<ParcelType>::NuclearEulerParcel
(
    const NuclearEulerParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh)
{}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "NuclearEulerParcelIO.C"

// ************************************************************************* //
