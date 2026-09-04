{
    // Generic loop over all delayed neutron precursor groups
    // (was 8 near-identical copy-pasted precEq1..precEq8 blocks)

    diffEff =(thermo2.nu()/Sc + phase2.nut()/Sct);//+(thermo1.nu()/Sc + phase1.turbulence().nut()/Sct)*(rho1*alpha1)/(rho1*alpha1+rho2*alpha2);

    forAll(prec, i)
    {
        volScalarField& preci = prec[i];

        fvScalarMatrix precEqni
        (
            fvm::ddt(preci)
          + fvm::div(phi2, preci)
          - fvm::laplacian(diffEff, preci)
          + fvm::Sp(precLambda[i], preci)
         ==
            (1.0/Keff)*(1.0/KeffControl)*NSFflu*precBeta[i]
          + fvOptions(preci)
        );

        precEqni.relax();
        fvOptions.constrain(precEqni);
        precEqni.solve();
    }

	for (int eI=0; eI<energyGroups ; ++eI)
	{
		delNeuSource[eI] = chiDel[eI] * precLambda[0] * prec[0] * 0;
		forAll(prec, j)
		{
			delNeuSource[eI] += chiDel[eI] * precLambda[j] * prec[j];
		}
    }
}
