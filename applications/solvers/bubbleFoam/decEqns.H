{
    //diffusion coefficient
    diffEff =(thermo2.nu()/Sc + phase2.nut()/Sct);//+(thermo1.nu()/Sc + phase1.turbulence().nut()/Sct)*(rho1*alpha1)/(rho1*alpha1+rho2*alpha2);

    // Generic loop over all decay heat precursor groups
    // (was 3 near-identical copy-pasted decEq1..decEq3 blocks)

    forAll(dec, i)
    {
        volScalarField& deci = dec[i];

        fvScalarMatrix decEqni
        (
            fvm::ddt(deci)
          + fvm::div(phi2, deci)
          - fvm::laplacian(diffEff, deci)
          + fvm::Sp(decLambda[i], deci)
         ==
            SPflu*decBeta[i]
          + fvOptions(deci)
        );

        decEqni.relax();
        fvOptions.constrain(decEqni);
        decEqni.solve();
    }

    decayHeatSource = dec[0]*decLambda[0]*0;
    forAll(dec, i)
    {
        decayHeatSource += decLambda[i]*dec[i];
    }
}
