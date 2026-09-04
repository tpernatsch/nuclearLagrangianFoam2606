/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "flowRateForce.H"
#include "fvMatrices.H"
#include "DimensionedField.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(flowRateForce, 0);

    addToRunTimeSelectionTable
    (
        option,
        flowRateForce,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::flowRateForce::writeProps
(
    const scalar gradP
) const
{
    // Only write on output time
    if (mesh_.time().writeTime())
    {
        IOdictionary propsDict
        (
            IOobject
            (
                name_ + "Properties",
                mesh_.time().timeName(),
                "uniform",
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );
        propsDict.add("gradient", gradP);
        propsDict.regIOobject::write();
    }
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::flowRateForce::initialise()
{
    const faceZone& fZone = mesh_.faceZones()[zoneID_];

    faceId_.setSize(fZone.size());
    facePatchId_.setSize(fZone.size());
    faceSign_.setSize(fZone.size());

    label count = 0;
    forAll(fZone, i)
    {
        label facei = fZone[i];
        label faceId = -1;
        label facePatchId = -1;
        if (mesh_.isInternalFace(facei))
        {
            faceId = facei;
            facePatchId = -1;
        }
        else
        {
            facePatchId = mesh_.boundaryMesh().whichPatch(facei);
            const polyPatch& pp = mesh_.boundaryMesh()[facePatchId];
            if (isA<coupledPolyPatch>(pp))
            {
                if (refCast<const coupledPolyPatch>(pp).owner())
                {
                    faceId = pp.whichFace(facei);
                }
                else
                {
                    faceId = -1;
                }
            }
            else if (!isA<emptyPolyPatch>(pp))
            {
                faceId = facei - pp.start();
            }
            else
            {
                faceId = -1;
                facePatchId = -1;
            }
        }

        if (faceId >= 0)
        {
            if (fZone.flipMap()[i])
            {
                faceSign_[count] = -1;
            }
            else
            {
                faceSign_[count] = 1;
            }
            faceId_[count] = faceId;
            facePatchId_[count] = facePatchId;
            count++;
        }
    }
    faceId_.setSize(count);
    facePatchId_.setSize(count);
    faceSign_.setSize(count);

    calculateTotalArea(faceZoneArea_);
}


void Foam::fv::flowRateForce::calculateTotalArea
(
    scalar& area
)
{
    area = 0;
    forAll(faceId_, i)
    {
        label facei = faceId_[i];
        if (facePatchId_[i] != -1)
        {
            label patchi = facePatchId_[i];
            area += mesh_.magSf().boundaryField()[patchi][facei];
        }
        else
        {
            area += mesh_.magSf()[facei];
        }
    }
    reduce(area, sumOp<scalar>());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::flowRateForce::flowRateForce
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(sourceName, modelType, dict, mesh),
    G_(readScalar(coeffs_.lookup("flowRate"))),
    gradP0_(0.0),
    dGradP_(0.0),
    flowVec_(coeffs_.lookup("direction")),
    flowDir_(flowVec_/mag(flowVec_)),
    relaxation_(coeffs_.lookupOrDefault<scalar>("relaxation", 1.0)),
    rAPtr_(nullptr),
    faceZoneName_(coeffs_.lookup("faceZone")),
    zoneID_(mesh_.faceZones().findZoneID(faceZoneName_)),
    faceId_(),
    facePatchId_(),
    faceSign_(),
    faceZoneArea_(0)
{
    coeffs_.lookup("fields") >> fieldNames_;

    if (fieldNames_.size() != 1)
    {
        FatalErrorInFunction
            << "settings are:" << fieldNames_ << exit(FatalError);
    }

    applied_.setSize(fieldNames_.size(), false);

    // Read the initial pressure gradient from file if it exists
    IFstream propsFile
    (
        mesh_.time().timePath()/"uniform"/(name_ + "Properties")
    );

    if (propsFile.good())
    {
        Info<< "    Reading pressure gradient from file" << endl;
        dictionary propsDict(dictionary::null, propsFile);
        propsDict.lookup("gradient") >> gradP0_;
    }

    Info<< "    Initial pressure gradient = " << gradP0_ << nl << endl;

    if (zoneID_ < 0)
    {
        FatalErrorInFunction
            << type() << " " << this->name() << ": "
            << "    Unknown face zone name: " << faceZoneName_
            << ". Valid face zones are: " << mesh_.faceZones().names()
            << nl << exit(FatalError);
    }

    initialise();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::fv::flowRateForce::magFlowRate() const
{
    const surfaceScalarField& phi =
        mesh_.lookupObject<surfaceScalarField>("phi");

    scalar totalphi = 0;
    forAll(faceId_, i)
    {
        label facei = faceId_[i];
        if (facePatchId_[i] != -1)
        {
            label patchi = facePatchId_[i];
            totalphi += phi.boundaryField()[patchi][facei]*faceSign_[i];
        }
        else
        {
            totalphi += phi[facei]*faceSign_[i];
        }
    }
    reduce(totalphi, sumOp<scalar>());

    return mag(totalphi);
}


void Foam::fv::flowRateForce::correct(volVectorField& U)
{
    const scalarField& rAU = rAPtr_();

    // Integrate flow variables over cell set
    scalar rAUave = 0.0;
    const scalarField& cv = mesh_.V();
    forAll(cells_, i)
    {
        label celli = cells_[i];
        scalar volCell = cv[celli];
        rAUave += rAU[celli]*volCell;
    }

    // Collect across all processors
    reduce(rAUave, sumOp<scalar>());

    // Volume averages
    rAUave /= V_;

    scalar G = this->magFlowRate();

    // Calculate the pressure gradient increment needed to adjust the
    // flow-rate to the desired value
    //dGradP_ = relaxation_*(mag(Ubar_) - magUbarAve)/rAUave;
    dGradP_ = relaxation_*(G_ - G)/(rAUave*faceZoneArea_);

    // Apply correction to velocity field
    forAll(cells_, i)
    {
        label celli = cells_[i];
        U[celli] += flowDir_*rAU[celli]*dGradP_;
    }

    scalar gradP = gradP0_ + dGradP_;

    Info<< "Pressure gradient source: uncorrected flow rate = " << G
        << ", pressure gradient = " << gradP << endl;

    writeProps(gradP);
}


void Foam::fv::flowRateForce::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    volVectorField::Internal Su
    (
        IOobject
        (
            name_ + fieldNames_[fieldi] + "Sup",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", eqn.dimensions()/dimVolume, Zero)
    );

    scalar gradP = gradP0_ + dGradP_;

    UIndirectList<vector>(Su, cells_) = flowDir_*gradP;

    eqn += Su;
}


void Foam::fv::flowRateForce::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    this->addSup(eqn, fieldi);
}


void Foam::fv::flowRateForce::constrain
(
    fvMatrix<vector>& eqn,
    const label
)
{
    if (!rAPtr_)
    {
        rAPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    name_ + ":rA",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                1.0/eqn.A()
            )
        );
    }
    else
    {
        rAPtr_() = 1.0/eqn.A();
    }

    gradP0_ += dGradP_;
    dGradP_ = 0.0;
}


bool Foam::fv::flowRateForce::read(const dictionary& dict)
{
    NotImplemented;

    return false;
}


// ************************************************************************* //
