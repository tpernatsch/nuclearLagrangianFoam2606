/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "NuclearParcel.H"
#include "physicoChemicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::NuclearParcel<ParcelType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    ParcelType::setCellValues(cloud, td);

    tetIndices tetIs = this->currentTetIndices();

    //td.Cpc() = cloud.constProps().Cp0();

    td.Tc() = td.TInterp().interpolate(this->coordinates(), tetIs);

    if (td.Tc() < cloud.constProps().TMin())
    {
        if (debug)
        {
            WarningInFunction
                << "Limiting observed temperature in cell " << this->cell()
                << " to " << cloud.constProps().TMin() <<  nl << endl;
        }

        td.Tc() = cloud.constProps().TMin();
    }
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::NuclearParcel<ParcelType>::cellValueSourceCorrection
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    const label celli = this->cell();
    const scalar massCell = this->massCell(td);

    td.Uc() += cloud.UTrans()[celli]/massCell;

    const scalar CpMean = cloud.constProps().Cp0();
    td.Tc() += cloud.hsTrans()[celli]/(CpMean*massCell);

    if (td.Tc() < cloud.constProps().TMin())
    {
        if (debug)
        {
            WarningInFunction
                << "Limiting observed temperature in cell " << celli
                << " to " << cloud.constProps().TMin() <<  nl << endl;
        }

        td.Tc() = cloud.constProps().TMin();
    }
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::NuclearParcel<ParcelType>::calcSurfaceValues
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar T,
    scalar& Ts,
    scalar& rhos,
    scalar& mus,
    scalar& Pr,
    scalar& kappas
) const
{
    // Surface temperature using two thirds rule
    Ts = (2.0*T + td.Tc())/3.0;

    if (Ts < cloud.constProps().TMin())
    {
        if (debug)
        {
            WarningInFunction
                << "Limiting parcel surface temperature to "
                << cloud.constProps().TMin() <<  nl << endl;
        }

        Ts = cloud.constProps().TMin();
    }

    // For incompressible: all properties are constant
    rhos = cloud.constProps().rho0();
    mus = cloud.constProps().mu0();
    kappas = cloud.constProps().kappa0();
    Pr = mus * cloud.constProps().Cp0() / kappas;
    Pr = max(ROOTVSMALL, Pr);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::NuclearParcel<ParcelType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    // Define local properties at beginning of time step
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar np0 = this->nParticle_;
    const scalar mass0 = this->mass();

    // Store T for consistent radiation source
    const scalar T0 = this->T_;


    // Calc surface values
    // ~~~~~~~~~~~~~~~~~~~
    scalar Ts, rhos, mus, Pr, kappas;
    calcSurfaceValues(cloud, td, this->T_, Ts, rhos, mus, Pr, kappas);

    // Reynolds number
    // Floored at SMALL (not just >0): a fully Stokes-equilibrated parcel
    // (slip velocity ~0) can underflow Re through denormals down to
    // exactly 0.0, which reliably SIGFPEs somewhere in the
    // heatTransferModel::htc()/Nu() call chain for AhmedYovanovich
    // despite its own Re<=0 entry guard (root cause not pinned down -
    // see project memory). Flooring here, at the source, means neither
    // Nu() nor anything else downstream ever sees the problematic value,
    // regardless of where exactly it was faulting.
    scalar Re = max(this->Re(rhos, this->U_, td.Uc(), this->d_, mus), SMALL);


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

    // Calculate new particle temperature
    this->T_ =
        this->calcHeatTransfer
        (
            cloud,
            td,
            dt,
            Re,
            Pr,
            kappas,
            NCpW,
            Sh,
            dhsTrans,
            Sph
        );


    // Motion
    // ~~~~~~
    
    //Kappa adjustment to take into account particle temperature
    //kappa *= T_**theta

    // Calculate new particle velocity
    this->U_ =
        this->calcVelocity(cloud, td, dt, Re, mus, mass0, Su, dUTrans, Spu); //add kappa for diffusion calculation



    //  Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
        cloud.hsCoeffTemp()[this->cell()] += np0*Sph*T_;

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


template<class ParcelType>
template<class TrackCloudType>
Foam::scalar Foam::NuclearParcel<ParcelType>::calcHeatTransfer
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const scalar Re,
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
        return T_;
    }

    const scalar d = this->d();
    const scalar rho = cloud.constProps().rho0();
    const scalar As = this->areaS(d);
    const scalar V = this->volume(d);
    const scalar m = rho*V;
    const scalar Cp = cloud.constProps().Cp0();

    // Calc heat transfer coefficient
    scalar htc = cloud.heatTransfer().htc(d, Re, Pr, kappa, NCpW);
    scalar decHeat = cloud.decayHeat().decayPower(d);  //Addition for deacying particles

    // Calculate the integration coefficients
    const scalar bcp = htc*As/(m*Cp);
    const scalar acp = bcp*td.Tc(); //+ decHeat/(m*Cp);  // decHeat/(m*Cp) is the addition for decaying particles

    scalar ancp = Sh;
    if (cloud.radiation())
    {
        const tetIndices tetIs = this->currentTetIndices();
        const scalar Gc = td.GInterp().interpolate(this->coordinates(), tetIs);
        const scalar sigma = physicoChemical::sigma.value();
        const scalar epsilon = cloud.constProps().epsilon0();

        ancp += As*epsilon*(Gc/4.0 - sigma*pow4(T_));
    }
    ancp += decHeat;  // Addition for decaying particles
    ancp /= m*Cp;

    // Integrate to find the new parcel temperature
    const scalar deltaT = cloud.TIntegrator().delta(T_, dt, acp + ancp, bcp);
    const scalar deltaTncp = ancp*dt;
    const scalar deltaTcp = deltaT - deltaTncp;

    // Calculate the new temperature and the enthalpy transfer terms
    scalar Tnew = T_ + deltaT;
    Tnew = clamp(Tnew, cloud.constProps().TMin(), cloud.constProps().TMax());

    dhsTrans -= m*Cp*deltaTcp;
    Sph = dt*m*Cp*bcp;
    

    return Tnew;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::NuclearParcel<ParcelType>::NuclearParcel
(
    const NuclearParcel<ParcelType>& p
)
:
    ParcelType(p),
    T_(p.T_),
    Cp_(p.Cp_)
{}


template<class ParcelType>
Foam::NuclearParcel<ParcelType>::NuclearParcel
(
    const NuclearParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    T_(p.T_),
    Cp_(p.Cp_)
{}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "NuclearParcelIO.C"

// ************************************************************************* //
