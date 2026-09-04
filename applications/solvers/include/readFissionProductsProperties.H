Info<< "\nReading fission products constants" << endl;

dictionary fissionProductsDict(nuclearPropertiesDict.subDict("fissionProducts"));

// Set pointer list size
nFPSpecies = readInt(fissionProductsDict.lookup("species"));
FPYield.setSize(nFPSpecies);
FPLambda.setSize(nFPSpecies);
FPDp.setSize(nFPSpecies);
FPRho.setSize(nFPSpecies);
FPPd.setSize(nFPSpecies);

// Scalar lists to lookup dictionary
List<scalar> FPYield_(fissionProductsDict.lookup("yield"));	//Fission yield
List<scalar> FPLambda_(fissionProductsDict.lookup("lambda"));	//Decay constant
List<scalar> FPDp_(fissionProductsDict.lookup("dp"));		//Particles diameter
List<scalar> FPRho_(fissionProductsDict.lookup("rho"));		//Particle density
List<scalar> FPPd_(fissionProductsDict.lookup("Pd"));		//Decay heat

// Check input data size
if( FPLambda_.size() == nFPSpecies && FPYield_.size() == nFPSpecies )
{
    Info<< "    Found " << nFPSpecies << " species of fission products" << endl;
}
else
{
    FatalError
        << "in " << fissionProductsDict.name() << " , invalid number of fission product species constants" << exit(FatalError);
}

// Assign dimensioned values to pointer list
forAll(FPYield_, i)
{
    FPYield.set(i, new dimensionedScalar("FPYield" + name(i + 1), dimless, FPYield_[i]));
    FPLambda.set(i, new dimensionedScalar("FPLambda" + name(i + 1), dimless/dimTime, FPLambda_[i]));
    FPDp.set(i, new dimensionedScalar("FPDp" + name(i + 1), dimLength, FPDp_[i]));
    FPRho.set(i, new dimensionedScalar("FPRho" + name(i + 1), dimDensity, FPRho_[i]));
    FPPd.set(i, new dimensionedScalar("FPPd" + name(i + 1), dimEnergy, FPPd_[i]));
}
