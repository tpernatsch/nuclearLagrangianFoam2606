{
    forAll(fp, i)
    {
        // Stokes-Einstein molecular diffusivity (via FPDp) + turbulent eddy diffusivity
        Dfp[i] =
            constant::physicoChemical::k*T
           /(3*constant::mathematical::pi*muc*FPDp[i])
          + turbulence->nut()/Sct;

        volScalarField& fpi = fp[i];
        volScalarField& Dfpi = Dfp[i];
        volScalarField& fpdi = fpd[i];

        fvScalarMatrix fpEqni
        (
            fvm::ddt(fpi)
          + fvm::div(phi, fpi)
          - fvm::laplacian(Dfpi, fpi)
          + fvm::Sp(FPLambda[i], fpi)
         ==
            FPYield[i]*Rf
          + fvOptions(fpi)
        );

        fpEqni.relax();
        fvOptions.constrain(fpEqni);
        fpEqni.solve();
        fvOptions.correct(fpi);

        #include "calcDeposition.H"
    }
}
