powDens = promptHeatSource + decayHeatSource;

scalar totalPower = fvc::domainIntegrate(powDens).value();

if(criticalityCalculation)
{
	scalar normalization = totalPower/nominalPower;

	for (int eI=0; eI<energyGroups ; ++eI)
	{
		flux[eI] = flux[eI]/normalization;

		if (neutronTransport)
		{
			secondFlux[eI] = secondFlux[eI]/normalization;
		}

		othGroFisNeu[eI]=othGroFisNeu[eI]/normalization;
		scaNeu[eI]=scaNeu[eI]/normalization;
		delNeuSource[eI]=delNeuSource[eI]/normalization;
	}

/*
	prec1 = prec1/normalization;
	prec2 = prec2/normalization;
	prec3 = prec3/normalization;
	prec4 = prec4/normalization;
	prec5 = prec5/normalization;
	prec6 = prec6/normalization;
	prec7 = prec7/normalization;
	prec8 = prec8/normalization;

	dec1 = dec1/normalization;
	dec2 = dec2/normalization;
	dec3 = dec3/normalization;
*/
	NSFflu=NSFflu/normalization;
	promptHeatSource=promptHeatSource/normalization;
	decayHeatSource=decayHeatSource/normalization;

	fissionRate=fissionRate/normalization;
}

Info << nl << "totalPower = " << fvc::domainIntegrate(powDens).value() << "W" << endl;
Info << nl << "fission power = " << fvc::domainIntegrate(promptHeatSource).value() << "W" << endl;
Info << nl << "void fraction = " << fvc::domainIntegrate(alpha1).value()/reactorVolume << endl;

//Integral inventory of Xe
Info << nl << "Xenon Inventory = " << fvc::domainIntegrate(xenonMass).value() << "kg" << endl;

NSFflu_old=NSFflu;

if (!criticalityCalculation)
{
	volPow=volPow/volPow*totalPower;
}
