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

Description
    Registers the secondary-Euler-phase dispersion models (see
    TwoPhaseDispersionRASModel and its derivatives) into
    NuclearEulerCloud::dispersion2()'s runtime-selection table - the SAME
    table as the primary phase's dispersionModel entry, since both are
    keyed by CloudType::kinematicCloudType (see makeDispersionModelType()
    in DispersionModel.H) and dispersion2()'s dict-key-rename trick (see
    NuclearEulerCloud.C) hands the "dispersionModel2" entry to the very
    same DispersionModel<kinematicCloudType>::New(). This only needs
    makeDispersionModelType() (as makeParcelTurbulenceDispersionModels.H
    does for the primary phase), never makeDispersionModel() again: the
    table itself is already defined once, in myIntermediate's own
    makeParcelDispersionModels(basicNuclearEulerCloud) call.

\*---------------------------------------------------------------------------*/

#ifndef makeNuclearEulerParcelTurbulenceDispersionModels_H
#define makeNuclearEulerParcelTurbulenceDispersionModels_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "TwoPhaseGradientDispersionRAS.H"
#include "TwoPhaseStochasticDispersionRAS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeNuclearEulerParcelTurbulenceDispersionModels(CloudType)           \
                                                                               \
    typedef Foam::CloudType::kinematicCloudType kinematicCloudType;          \
    defineNamedTemplateTypeNameAndDebug                                      \
    (                                                                        \
        Foam::TwoPhaseDispersionRASModel<kinematicCloudType>,                \
        0                                                                    \
    );                                                                       \
                                                                               \
    makeDispersionModelType(TwoPhaseGradientDispersionRAS, CloudType);       \
    makeDispersionModelType(TwoPhaseStochasticDispersionRAS, CloudType);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
