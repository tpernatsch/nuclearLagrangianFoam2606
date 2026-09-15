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

Class
    Foam::TwoPhaseDispersionRASModel

Description
    Base class for the secondary-Euler-phase dispersion models plugged into
    NuclearEulerCloud::dispersion2() (see nuclearEulerCloud.H and
    NuclearEulerParcel::calcDispersion()).

    Mirrors the role DispersionRASModel plays for the primary phase, but
    the mechanics of getting at the secondary phase's turbulence fields are
    different, and simpler: DispersionRASModel::kModel()/epsilonModel()
    must locate the primary phase's turbulence model in the mesh object
    registry and return it as a tmp<> (it may or may not already be a
    registered field), because plain KinematicCloud has no turbulence
    awareness of its own. NuclearEulerCloud, by contrast, already stores
    the secondary phase's k2()/epsilon2() as persistent fields, refreshed
    once per PIMPLE iteration by nuclearEulerCloudSource - so k2Field()/
    epsilon2Field() below just return references to those, no tmp<>/
    ownership bookkeeping required.

    The one piece of unavoidable bookkeeping is the same static-owner-type
    gotcha documented throughout NuclearEulerCloud/NuclearEulerParcel:
    DispersionModel<CloudType>::owner() is declared to return CloudType&,
    and CloudType here is fixed, by the makeDispersionModelType() macro, to
    CloudType::kinematicCloudType - i.e. the *KinematicCloud* layer of the
    template stack, which knows nothing about NuclearEulerCloud's own
    k2()/epsilon2() accessors. The object actually pointed to at runtime is
    always a basicNuclearEulerCloud (this dispersion model can only ever be
    selected for a NuclearEulerCloud's dispersion2() slot - see
    makeNuclearEulerParcelTurbulenceDispersionModels.H), so an explicit
    static_cast down to that concrete, non-templated type is safe and is
    exactly what the two-phase ParticleForces in myIntermediate already do
    to reach the same fields from their own owner().

SourceFiles
    TwoPhaseDispersionRASModel.C

\*---------------------------------------------------------------------------*/

#ifndef TwoPhaseDispersionRASModel_H
#define TwoPhaseDispersionRASModel_H

#include "DispersionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class TwoPhaseDispersionRASModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class TwoPhaseDispersionRASModel
:
    public DispersionModel<CloudType>
{
protected:

    // Protected Functions

        //- Return const access to the secondary phase's turbulence k
        const volScalarField& k2Field() const;

        //- Return const access to the secondary phase's turbulence
        //- epsilon
        const volScalarField& epsilon2Field() const;


public:

    //- Runtime type information
    TypeName("twoPhaseDispersionRASModel");


    // Constructors

        //- Construct from components
        TwoPhaseDispersionRASModel(const dictionary& dict, CloudType& owner);

        //- Construct copy
        TwoPhaseDispersionRASModel
        (
            const TwoPhaseDispersionRASModel<CloudType>& dm
        );

        //- Construct and return a clone
        virtual autoPtr<DispersionModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~TwoPhaseDispersionRASModel();


    // Member Functions

        //- Update (disperse particles)
        virtual vector update
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
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TwoPhaseDispersionRASModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
