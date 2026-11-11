/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

Application
    particleConcentrationFoam

Description
    Eulerian (concentration) counterpart of the Lagrangian particle runs of
    the heated pool, on a frozen flow.

    For every particle diameter d_i it solves, on the converged U, phi and T
    read from the start time,

        ddt(N_i) + div(phi, N_i) = laplacian(D_i, N_i) + S

    with the Brownian diffusivity of the Lagrangian BrownianMotion force,

        D_i = k_B T C_c / (3 pi mu d_i),  C_c = 1 + a (1.257 + 0.4 exp(-1.1/a)),
        a = 2 lambda / d_i,

    T the local carrier temperature and mu either constant (viscosity
    constant: as the parcels of runs without muFromTransport) or rho nu(T) of
    the case's own transport model (viscosity transport: as muFromTransport).

    S is a uniform volumetric source: the Eulerian form of the Lagrangian
    replacement of deposited parcels. The sticky patches take N = 0 (perfect
    sinks, set in 0/N). With ddtSchemes steadyState every time step is one
    outer iteration of the steady problem (under-relaxed with the
    relaxationFactors of N, deferred-correction schemes such as linearUpwind
    keep the implicit part an M-matrix); with Euler it marches in time. It
    stops when both the largest relative change of N in a step and the
    initial residual are below the tolerance.

    At every write it reports, per diameter and patch, the diffusive flux
    into the patch J = -D dN/dn, the mean concentration N_mean, the
    deposition velocity v_d = (integral J dA / A) / N_mean and the balance
    (sum of the patch fluxes)/(S V), and it writes the flux into the patches
    as the boundary values of the field J_<label> (deposition profiles).

    Input: constant/particleConcentrationProperties, e.g.

        diameters   (1e-9 2e-9 5e-9 1e-8 2e-8);   // [m]
        labels      (1nm 2nm 5nm 10nm 20nm);       // field names N_<label>
        viscosity   constant;       // constant | transport
        mu          0.00450913;     // [Pa s], viscosity constant
        rho         4364.9;         // [kg/m3], viscosity transport: mu = rho nu
        lambda      1e-12;          // [m], Cunningham (as BrownianMotion)
        source      1;              // [1/m3/s]
        tolerance   1e-7;           // steady when max |dN|/max N < tolerance
        phiScale    1;              // verification only: 0 = pure diffusion

    and a template 0/N (the boundary conditions of every N_<label>).

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "physicoChemicalConstants.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady particle concentration with Brownian diffusion on a frozen"
        " flow (Eulerian counterpart of the Lagrangian heated-pool runs)"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOdictionary props
    (
        IOobject
        (
            "particleConcentrationProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const scalarList diameters(props.get<scalarList>("diameters"));
    const wordList labels(props.get<wordList>("labels"));
    if (labels.size() != diameters.size())
    {
        FatalErrorInFunction
            << "diameters and labels differ in length" << exit(FatalError);
    }
    const word viscosity(props.get<word>("viscosity"));
    const scalar lambda(props.getOrDefault<scalar>("lambda", 1e-12));
    const dimensionedScalar S
    (
        "S", dimless/dimVolume/dimTime, props.getOrDefault<scalar>("source", 1)
    );
    const scalar tolerance(props.getOrDefault<scalar>("tolerance", 1e-7));

    Info<< "Reading U, phi and T" << endl;
    volVectorField U
    (
        IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
        mesh
    );
    surfaceScalarField phi
    (
        IOobject("phi", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
        mesh
    );
    volScalarField T
    (
        IOobject("T", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
        mesh
    );

    // phiScale: verification only (0 = pure diffusion, which has exact solutions)
    const scalar phiScale(props.getOrDefault<scalar>("phiScale", 1));
    if (phiScale != 1)
    {
        Info<< "phi scaled by " << phiScale << " (verification)" << endl;
        phi *= phiScale;
    }

    // Carrier dynamic viscosity: constant, or rho nu(T) of the transport model
    volScalarField mu
    (
        IOobject("mu", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar(dimDynamicViscosity, Zero)
    );
    autoPtr<singlePhaseTransportModel> transport;
    if (viscosity == "constant")
    {
        mu = dimensionedScalar("mu", dimDynamicViscosity, props.get<scalar>("mu"));
    }
    else if (viscosity == "transport")
    {
        transport.reset(new singlePhaseTransportModel(U, phi));
        transport->correct();
        const dimensionedScalar rho("rho", dimDensity, props.get<scalar>("rho"));
        mu = rho*transport->nu();
    }
    else
    {
        FatalErrorInFunction
            << "viscosity must be constant or transport, not " << viscosity
            << exit(FatalError);
    }
    Info<< "mu: " << gMin(mu.primitiveField()) << " - " << gMax(mu.primitiveField())
        << " Pa s, T: " << gMin(T.primitiveField()) << " - "
        << gMax(T.primitiveField()) << " K" << nl << endl;

    const dimensionedScalar kB
    (
        "kB", dimEnergy/dimTemperature, constant::physicoChemical::k.value()
    );

    // Template with the boundary conditions of every concentration field
    const volScalarField Ntemplate
    (
        IOobject("N", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
        mesh
    );

    PtrList<volScalarField> N(diameters.size());
    PtrList<volScalarField> D(diameters.size());
    forAll(diameters, i)
    {
        const scalar d = diameters[i];
        const scalar a = 2*lambda/d;
        const scalar cc = 1 + a*(1.257 + 0.4*Foam::exp(-1.1/a));

        N.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "N_" + labels[i],
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                Ntemplate
            )
        );
        D.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "D_" + labels[i],
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                kB*T*cc/(constant::mathematical::pi*3*mu*dimensionedScalar(dimLength, d))
            )
        );
        Info<< "N_" << labels[i] << ": d " << d << " m, C_c " << cc
            << ", D " << gMin(D[i].primitiveField()) << " - "
            << gMax(D[i].primitiveField()) << " m2/s" << endl;
    }

    // Sticky patches: those where every N has a fixed value (zero)
    const fvPatchList& patches = mesh.boundary();
    DynamicList<label> sinks;
    forAll(patches, patchi)
    {
        if (Ntemplate.boundaryField()[patchi].fixesValue())
        {
            sinks.append(patchi);
        }
    }

    const scalar V = gSum(mesh.V());

    autoPtr<OFstream> logPtr;
    if (Pstream::master())
    {
        mkDir(runTime.globalPath()/"postProcessing");
        logPtr.reset
        (
            new OFstream(runTime.globalPath()/"postProcessing"/"depositionVelocity.dat")
        );
        logPtr() << "# time size patch area_m2 flux_per_s Nmean_per_m3 v_d_m_s"
                 << " balance" << endl;
    }

    Info<< nl << "Marching to the steady state" << nl << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << endl;

        scalar maxChange = 0;
        scalar maxResidual = 0;

        forAll(N, i)
        {
            volScalarField& Ni = N[i];
            const scalarField N0(Ni.primitiveField());

            fvScalarMatrix NEqn
            (
                fvm::ddt(Ni)
              + fvm::div(phi, Ni, "div(phi,N)")
              - fvm::laplacian(D[i], Ni, "laplacian(D,N)")
             ==
                S
            );
            NEqn.relax();
            maxResidual = max(maxResidual, NEqn.solve("N").initialResidual());

            const scalar dN = gMax(mag(Ni.primitiveField() - N0)());
            const scalar Nmax = max(gMax(Ni.primitiveField()), VSMALL);
            maxChange = max(maxChange, dN/Nmax);
        }

        Info<< "    max relative change " << maxChange
            << ", max initial residual " << maxResidual << endl;

        const bool steady = maxChange < tolerance && maxResidual < tolerance;

        if (runTime.writeTime() || steady)
        {
            forAll(N, i)
            {
                const volScalarField& Ni = N[i];
                const scalar Nmean = gSum(Ni.primitiveField()*mesh.V())/V;

                volScalarField J
                (
                    IOobject
                    (
                        "J_" + labels[i],
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedScalar(dimless/dimArea/dimTime, Zero),
                    fvPatchFieldBase::calculatedType()
                );

                scalar total = 0;
                Info<< "    " << labels[i] << ": N_mean " << Nmean;
                forAll(sinks, k)
                {
                    const label patchi = sinks[k];
                    const scalarField flux
                    (
                       -D[i].boundaryField()[patchi]
                       *Ni.boundaryField()[patchi].snGrad()
                    );
                    J.boundaryFieldRef()[patchi] = flux;

                    const scalarField& magSf = mesh.magSf().boundaryField()[patchi];
                    const scalar A = gSum(magSf);
                    const scalar F = gSum(flux*magSf);
                    total += F;
                    const scalar vd = (A > 0 && Nmean > 0) ? F/A/Nmean : 0;

                    Info<< ", " << patches[patchi].name() << " v_d " << vd;
                    if (Pstream::master())
                    {
                        logPtr() << runTime.timeName() << ' ' << labels[i] << ' '
                            << patches[patchi].name() << ' ' << A << ' ' << F
                            << ' ' << Nmean << ' ' << vd << ' '
                            << total/(S.value()*V) << endl;
                    }
                }
                Info<< " m/s, balance " << total/(S.value()*V) << endl;
                J.write();
            }
            runTime.writeNow();
        }

        if (steady)
        {
            Info<< nl << "Steady: max relative change " << maxChange
                << " < " << tolerance << nl << endl;
            break;
        }

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
