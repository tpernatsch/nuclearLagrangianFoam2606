/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoPhysicsTypes.H"

#include "rhoConst.H"
#include "perfectFluid.H"
#include "Boussinesq.H"

#include "pureMixture.H"
#include "multiComponentMixture.H"
#include "reactingMixture.H"
#include "SpecieMixture.H"

#include "rhoThermo.H"
#include "rhoReactionThermo.H"
#include "heRhoThermo.H"

#include "solidThermo.H"
#include "heSolidThermo.H"
#include "solidThermoPhysicsTypes.H"

#include "gasDiffusionSimple.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef
        constTransport
        <
            species::thermo
            <
                hConstThermo
                <
                    rhoConst<specie>
                >,
                sensibleEnthalpy
            >
        > constRhoHThermoPhysics;


    typedef
        constTransport
        <
            species::thermo
            <
                hConstThermo
                <
                    Boussinesq<specie>
                >,
                sensibleEnthalpy
            >
        > BoussinesqFluidEThermoPhysics;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    using namespace meltingEvaporationModels;

    //NOTE: First thermo (from) and second otherThermo (to)

        // From multi-component liquid (rhoConst) to a multi-component gas incomp phase
        makeSpecieInterfaceSpecieMixtures
        (
            gasDiffusionSimple,
            heRhoThermo,
            rhoReactionThermo,
            multiComponentMixture,
            constRhoHThermoPhysics,
            heRhoThermo,
            rhoReactionThermo,
            multiComponentMixture,
            constIncompressibleGasHThermoPhysics
        );

        // From multi-component liquid (rhoConst) to a multi-component liquid (rhoConst)
        makeSpecieInterfaceSpecieMixtures
        (
            gasDiffusionSimple,
            heRhoThermo,
            rhoReactionThermo,
            multiComponentMixture,
            constRhoHThermoPhysics,
            heRhoThermo,
            rhoReactionThermo,
            multiComponentMixture,
            constRhoHThermoPhysics
        );

}


// ************************************************************************* //
