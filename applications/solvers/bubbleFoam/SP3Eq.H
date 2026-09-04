{
	#include "updateXS.H"

	for (int eI=0; eI<energyGroups ; ++eI)
	{
		diffCoeffAlbedo=1.0/(3.0*TR[eI]);
		fluxAlbedo=secondFlux[eI];
		flux[eI].correctBoundaryConditions();

		fvScalarMatrix NEqn
		(
	        (1-criticalityCalculation)*fvm::ddt(IV[eI],flux[eI])
	  	  - fvm::laplacian(1.0/(3.0*TR[eI]),flux[eI])
		  - 1.0/Keff/KeffControl*fvm::Sp(F[eI],flux[eI])
	      + 2.0/Keff/KeffControl*fvc::Sp(F[eI],secondFlux[eI])
	      + fvm::Sp(M[eI],flux[eI])
	      - 2.0*fvc::Sp(M[eI],secondFlux[eI])
	      - delNeuSource[eI]
	      - scaNeu[eI]
	      - 1.0/Keff/KeffControl*othGroFisNeu[eI]
	      - (1-criticalityCalculation)*2.0*fvc::ddt(IV[eI],secondFlux[eI])
		);

		NEqn.relax();
		fvOptions.constrain(NEqn);
		NEqn.solve();
	}

	for (int eI=0; eI<energyGroups ; ++eI)
	{
        diffCoeffAlbedo=9.0/35.0/TR[eI];
        //diffCoeffAlbedo=9.0/35.0/(M[eI]+S[eI][eI]);
        //diffCoeffAlbedo=3.0/7.0/(M[eI]+S[eI][eI]);
    	fluxAlbedo=flux[eI];
		secondFlux[eI].correctBoundaryConditions();

    	fvScalarMatrix secondNEqn
    	(
            (1-criticalityCalculation)*9.0/5.0*fvm::ddt(IV[eI],secondFlux[eI])
      	  - 9.0/35.0*fvm::laplacian(1.0/TR[eI],secondFlux[eI])
      	  + fvm::Sp(TR[eI],secondFlux[eI])
      	  - 2.0/5.0*fvc::Sp(M[eI],flux[eI])
      	  + 4.0/5.0*fvm::Sp(M[eI],secondFlux[eI])
      	  + 2.0/5.0*1.0/Keff/KeffControl*fvc::Sp(F[eI],flux[eI])
          - 4.0/5.0*1.0/Keff/KeffControl*fvm::Sp(F[eI],secondFlux[eI])
      	  + 2.0/5.0*1.0/Keff/KeffControl*othGroFisNeu[eI]
      	  + 2.0/5.0*delNeuSource[eI]
      	  + 2.0/5.0*scaNeu[eI]
      	  - (1-criticalityCalculation)*2.0/5.0*fvc::ddt(IV[eI],flux[eI])

    	);

    	secondNEqn.relax();
    	fvOptions.constrain(secondNEqn);
    	secondNEqn.solve();
	}

	//azzera termini di sorgente
    for (int eI=0.0; eI<energyGroups ; ++eI)
    {
        scaNeu[eI]*=0.0;
        othGroFisNeu[eI]*=0.0;
    }

	NSFflu *= 0.0;
    SPflu *= 0.0;

	//aggiorna termini di sorgente
    for (int eI=0; eI<energyGroups ; ++eI)
    {
        NSFflu = NSFflu + NSF[eI] * (flux[eI] - 2.0*secondFlux[eI]);
        SPflu  = SPflu  + SP[eI] * (flux[eI] - 2.0*secondFlux[eI]);
        // delNeuSource is (re)computed once per timestep, at the end of
        // precEqns.H after the precursor fields are solved - no need to
        // duplicate that computation here (it was previously copy-pasted
        // in three places using the same stale prec fields)

		for (int eJ=0; eJ<energyGroups ; ++eJ)
        {
            if (eI!=eJ)
            {
                scaNeu[eI] = scaNeu[eI] + S[eJ][eI] * (flux[eJ] - 2.0*secondFlux[eJ]);
                othGroFisNeu[eI] = othGroFisNeu[eI] + NSF[eJ] * (1.0 - precBetaTot) * (flux[eJ] - 2.0*secondFlux[eJ])*chiPro[eI];
            }
        }
    }

    promptHeatSource = SPflu;

		// --- Xenon production/removal
			fissionRate = NSFflu / nubar;

			xeCapture *= 0.0;
			for (int eI=0; eI<energyGroups; ++eI)
			{
				xeCapture = xeCapture + alfaXE[eI] * flux[eI];
			}
}
