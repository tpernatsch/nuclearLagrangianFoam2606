{
    volScalarField Ddec("Ddec", turbulence->nu()/Sc + turbulence->nut()/Sct);

    forAll(dec, i)
    {
        volScalarField& deci = dec[i];

        fvScalarMatrix decEqni
        (
            fvm::ddt(deci)
          + fvm::div(phi, deci)
          - fvm::laplacian(Ddec, deci)
          + fvm::Sp(decLambda[i], deci)
         ==
            decBeta[i]*EfRf
          + fvOptions(deci)
        );

        decEqni.relax();
	fvOptions.constrain(decEqni);
        decEqni.solve();
    }
}
