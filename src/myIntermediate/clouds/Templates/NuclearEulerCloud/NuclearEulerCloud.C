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

#include "NuclearEulerCloud.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NuclearEulerCloud<CloudType>::NuclearEulerCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& muc,
    const volScalarField& T,
    const dimensionedVector& g,
    const volScalarField& alpha1,
    const volScalarField& k1,
    const volScalarField& epsilon1,
    const bool k1IsRAS,
    const eulerCarrierPhaseFields& phase2,
    bool readFields
)
:
    CloudType(cloudName, rho, U, muc, T, g, readFields),
    nuclearEulerCloud(),
    cloudCopyPtr_(nullptr),
    alpha1_(alpha1),
    alpha2_(phase2.alpha),
    U2_(phase2.U),
    rho2_(phase2.rho),
    mu2_(phase2.mu),
    T2_(phase2.T),
    k1_(k1),
    epsilon1_(epsilon1),
    k2_(phase2.k),
    epsilon2_(phase2.epsilon),
    d2_(phase2.d),
    k1IsRAS_(k1IsRAS),
    k2IsRAS_(phase2.isRAS),
    dispersion2_(nullptr),
    chiSwitchModel_(nullptr)
{
    // dispersionModel2 selects the secondary-phase dispersion model,
    // reusing DispersionModel::New()'s hardcoded "dispersionModel" key
    // via a renamed copy of subModelProperties(); defaults to "none".
    dictionary dispersion2Dict(this->subModelProperties());

    if (dispersion2Dict.found("dispersionModel2"))
    {
        dispersion2Dict.set
        (
            "dispersionModel",
            word(dispersion2Dict.lookup("dispersionModel2"))
        );
    }
    else
    {
        dispersion2Dict.set("dispersionModel", word("none"));
    }

    dispersion2_.reset
    (
        DispersionModel<typename CloudType::kinematicCloudType>::New
        (
            dispersion2Dict,
            *this
        ).ptr()
    );

    chiSwitchModel_.reset
    (
        ChiSwitchModel<nuclearEulerCloudType>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
}


template<class CloudType>
Foam::NuclearEulerCloud<CloudType>::NuclearEulerCloud
(
    NuclearEulerCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    nuclearEulerCloud(),
    cloudCopyPtr_(nullptr),
    alpha1_(c.alpha1_),
    alpha2_(c.alpha2_),
    U2_(c.U2_),
    rho2_(c.rho2_),
    mu2_(c.mu2_),
    T2_(c.T2_),
    k1_(c.k1_),
    epsilon1_(c.epsilon1_),
    k2_(c.k2_),
    epsilon2_(c.epsilon2_),
    d2_(c.d2_),
    k1IsRAS_(c.k1IsRAS_),
    k2IsRAS_(c.k2IsRAS_),
    dispersion2_(c.dispersion2_->clone()),
    chiSwitchModel_(c.chiSwitchModel_->clone())
{}


template<class CloudType>
Foam::NuclearEulerCloud<CloudType>::NuclearEulerCloud
(
    const fvMesh& mesh,
    const word& name,
    const NuclearEulerCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    nuclearEulerCloud(),
    cloudCopyPtr_(nullptr),
    alpha1_(c.alpha1_),
    alpha2_(c.alpha2_),
    U2_(c.U2_),
    rho2_(c.rho2_),
    mu2_(c.mu2_),
    T2_(c.T2_),
    k1_(c.k1_),
    epsilon1_(c.epsilon1_),
    k2_(c.k2_),
    epsilon2_(c.epsilon2_),
    d2_(c.d2_),
    k1IsRAS_(c.k1IsRAS_),
    k2IsRAS_(c.k2IsRAS_),
    dispersion2_(nullptr),
    chiSwitchModel_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::NuclearEulerCloud<CloudType>::preEvolve
(
    const typename parcelType::trackingData& td
)
{
    CloudType::preEvolve(td);

    // Cache secondary-phase dispersion fields, RAS phases only
    if (k2IsRAS_)
    {
        dispersion2_->cacheFields(true);
    }
}


template<class CloudType>
void Foam::NuclearEulerCloud<CloudType>::postEvolve
(
    const typename parcelType::trackingData& td
)
{
    CloudType::postEvolve(td);

    if (k2IsRAS_)
    {
        dispersion2_->cacheFields(false);
    }
}


template<class CloudType>
void Foam::NuclearEulerCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::trackingData td(*this);

        this->solve(*this, td);
    }
}


template<class CloudType>
void Foam::NuclearEulerCloud<CloudType>::autoMap(const mapPolyMesh& mapper)
{
    Cloud<parcelType>::autoMap(mapper);

    this->updateMesh();
}


template<class CloudType>
void Foam::NuclearEulerCloud<CloudType>::info()
{
    CloudType::info();

    Log_<< "    Secondary-phase volume fraction min/max = "
        << gMin(alpha2_.primitiveField())
        << ", " << gMax(alpha2_.primitiveField()) << nl
        << "    Secondary-phase temperature min/max     = "
        << gMin(T2_.primitiveField())
        << ", " << gMax(T2_.primitiveField()) << endl;
}


// ************************************************************************* //
