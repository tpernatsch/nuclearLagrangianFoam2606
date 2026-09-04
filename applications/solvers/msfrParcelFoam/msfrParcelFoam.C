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

Application
    msfrParcelFoam

Description
    Transient solver for buoyant, turbulent flow of incompressible fluids.
    Expanded to include neutronics and nuclear source terms for MSFR applications.
    Complete coupling of flow, temperature, and neutronics.
    Addition of lagranian particle tracking for noble metals.

    Uses the Boussinesq approximation:
    \f[
        rho_{k} = 1 - beta(T - T_{ref})
    \f]

    where:
        \f$ rho_{k} \f$ = the effective (driving) kinematic density
        beta = thermal expansion coefficient [1/K]
        T = temperature [K]
        \f$ T_{ref} \f$ = reference temperature [K]

    Valid when:
    \f[
        \frac{beta(T - T_{ref})}{rho_{ref}} << 1
    \f]

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "surfaceFilmModel.H"
#include "basicNuclearCloud.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "radiationModel.H"
#include "localEulerDdtScheme.H"
#include "physicoChemicalConstants.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

	argList::addNote
    (
        "Transient solver for incompressible, turbulent flow"
        " with nuclear particle clouds"
        " and surface film modelling. For the thermic part the" 
        " Boussinesq approximation is used."
        " The particles are modelled as nuclear particles,"
        " i.e. they undergo decay, and emit decay heat."
        " The neutronics are modelled using a multi-group diffusion approach,"
        " and the power distribution is updated at each time step."
    );

    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createNuclearFields.H"
    #include "createRegionControls.H"
    #include "createUfIfPresent.H"

    turbulence->validate();

    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        FP.storeGlobalPositions();
        FP.evolve();
        
        surfaceFilm.evolve();

        if (solvePrimaryRegion)
        {
        // --- Pressure-velocity PIMPLE corrector loop
            while (pimple.loop())
            {
                #include "UEqn.H"
                vort = fvc::curl(U);

                #include "TEqn.H"

                // --- Pressure corrector loop
                while (pimple.correct())
                {
                    #include "pEqn.H"
                }

                if (pimple.turbCorr())
                {
                    laminarTransport.correct();
                    turbulence->correct();
                }

                #include "updateCrossSections.H"

                // --- Neutronics NPIMPLE loop
                while (npimple.loop())
                {
                    #include "fluxEqns.H"
                    #include "precEqns.H"
                    #include "decEqns.H"

                    #include "updateNeutronSource.H"
                    #include "updateFissionRate.H"
                }

                #include "updatePowerSource.H"

                #include "fpEqns.H"
            }
        }
        runTime.write();

        Info<< "Q (MW)      : " << max(Q).value()/1E+06 << endl;

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
