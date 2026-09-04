{
    volScalarField Dprec("Dprec", turbulence->nu()/Sc + turbulence->nut()/Sct);

    forAll(prec, i)
    {
        volScalarField& preci = prec[i];

        fvScalarMatrix precEqni
        (
            fvm::ddt(preci)
          + fvm::div(phi, preci)
          - fvm::laplacian(Dprec, preci)
          + fvm::Sp(precLambda[i], preci)
         ==
            precBeta[i]*NuRf
          + fvOptions(preci)
        );

        precEqni.relax();
	fvOptions.constrain(precEqni);
        precEqni.solve();
    }
}
