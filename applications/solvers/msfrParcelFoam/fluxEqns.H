{
    forAll(flux, i)
    {
        volScalarField& fluxi = flux[i];

        diffCoeffAlbedo = D[i];

        fluxi.correctBoundaryConditions();

        fvScalarMatrix fluxEqni
        (
            fvm::ddt(invVel[i], fluxi)
          - fvm::laplacian(D[i], fluxi)
          - fvm::Sp(M[i], fluxi)
         ==
            neutronSource[i]
        );

        fluxEqni.relax();
        fvOptions.constrain(fluxEqni);
        fluxEqni.solve();
    }

    // Update total flux
    fluxTot *= 0.0;
    forAll(flux, i)
    {
        fluxTot += flux[i];
    }
}
