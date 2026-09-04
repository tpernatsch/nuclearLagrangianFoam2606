{
	#include "updateXS.H"

	for (int eI=0; eI<energyGroups ; ++eI)
	{
		diffCoeffAlbedo=1.0/(3.0*TR[eI]);
		fluxAlbedo=nullFlux*0;
		flux[eI].correctBoundaryConditions();

		fvScalarMatrix NEqn
		(
			(1-criticalityCalculation)*fvm::ddt(IV[eI],flux[eI])
		  - fvm::laplacian(1.0/(3.0*TR[eI]),flux[eI])
		  + fvm::Sp(M[eI],flux[eI])
			//ADDING XENON EFFECT
//			+ fvm::Sp(XE[eI],flux[eI])
			/////////////////////
		  - (1.0/Keff/KeffControl)*fvm::Sp(F[eI],flux[eI])
		  - delNeuSource[eI]
		  - scaNeu[eI]
		  - (1.0/Keff/KeffControl)*othGroFisNeu[eI]
	  	);

		NEqn.relax();
		fvOptions.constrain(NEqn);
		NEqn.solve();
	}

	//azzera termini di sorgente
	for (int eI=0; eI<energyGroups ; ++eI)
	{
		scaNeu[eI]*=0.0;
  		othGroFisNeu[eI]*=0.0;
	}

	NSFflu *= 0.0;
    SPflu *= 0.0;

	//calcola nuovi termini di sorgente
	for (int eI=0; eI<energyGroups ; ++eI)
	{
		NSFflu = NSFflu + NSF[eI] * flux[eI];
		SPflu  = SPflu  + SP[eI] * flux[eI];

		for (int eJ=0; eJ<energyGroups ; ++eJ)
		{
			if (eI!=eJ)
			{
				scaNeu[eI] = scaNeu[eI] + S[eJ][eI] * flux[eJ];
				othGroFisNeu[eI] = othGroFisNeu[eI] + (1-precBetaTot) * NSF[eJ] * flux[eJ] * chiPro[eI];
			}
		}
	}

	promptHeatSource = (1-decBetaTot) * SPflu;

	// --- Xenon production/removal
		fissionRate = NSFflu / nubar;

		xeCapture *= 0.0;
		for (int eI=0; eI<energyGroups; ++eI)
		{
			xeCapture = xeCapture + alfaXE[eI] * flux[eI];
		}
}
