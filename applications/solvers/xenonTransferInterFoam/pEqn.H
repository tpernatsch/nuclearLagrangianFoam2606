{
    if (correctPhi)
    {
        rAU.ref() = 1.0/UEqn.A();
    }
    else
    {
        rAU = 1.0/UEqn.A();
    }

    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));

    volVectorField HbyA("HbyA", U);
    HbyA = rAU()*UEqn.H();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh.Sf())
      + fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi)
    );

    if (p_rgh.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p_rgh);
        fvc::makeAbsolute(phiHbyA, U);
    }

    surfaceScalarField phig
    (
        (
            stf
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the fixedFluxPressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf);

    while (pimple.correctNonOrthogonal())
    {

        fvScalarMatrix p_rghEqn
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
        );

        if (fluid.includeVolChange())
        {
            p_rghEqn += fluid.volTransfer(p_rgh);
        }

        p_rghEqn.setReference(pRefCell, pRefValue);

        p_rghEqn.solve(p_rgh.select(pimple.finalInnerIter()));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + p_rghEqn.flux();

            p_rgh.relax();

            surfaceScalarField fU((phig + p_rghEqn.flux())/rAUf);

            // Per drop, zero net (CSF - grad p_rgh) momentum (dropNetForce.H);
            // neglects the gas pressure force (buoyancy ~1e-6 g, form drag
            // < 1 % g here). Fluxes unchanged
            if (dropForceBalance != "none")
            {
                const dropNetForce dnf
                (
                    mesh.lookupObject<volScalarField>("alpha." + dropForceBalance),
                    rho
                );
                if (dnf.nDrops() > 0)
                {
                    // Volume weights: the conservative momentum balance
                    // (rho*rAU*V/dt weights left ~7 % of g, kin_fb 2026-11-23)
                    const vectorField T(dnf.nDrops(), Zero);
                    vectorField F0;
                    fU += dnf.balance
                    (
                        fU + ghf*fvc::snGrad(rho)*mesh.magSf(), mesh.V(), T, F0
                    );
                    if (pimple.finalInnerIter())
                    {
                        forAll(F0, j)
                        {
                            if (dnf.vol()[j] <= VSMALL) continue;
                            const scalar m = max(dnf.mass()[j], VSMALL);
                            Info<< "dropForceBalance drop " << j
                                << " net (csf - grad p_rgh)/m = " << F0[j]/m
                                << " m/s2" << endl;
                        }
                    }
                }
            }

            U = HbyA + rAU()*fvc::reconstruct(fU);

            U.correctBoundaryConditions();
            fvOptions.correct(U);
            K = 0.5*magSqr(U);
        }
    }

    #include "continuityErrs.H"

    // Correct Uf if the mesh is moving
    fvc::correctUf(Uf, U, phi);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }

    if (!correctPhi)
    {
        rAU.clear();
    }
}
