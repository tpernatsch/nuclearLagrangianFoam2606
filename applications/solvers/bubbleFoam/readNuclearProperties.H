Info << "\nReading neutron transport constants" << endl;

dictionary neutronTransportDict(nuclearProperties.subDict("neutronTransport"));

// Number of energy groups is nuclear data, not a run-control switch, so it
// lives in nuclearProperties (matches msfrParcelFoam's readNuclearProperties.H)
label energyGroups(readInt(neutronTransportDict.lookup("energyGroups")));

// DECLARATIONS

// Removal [1/m]
PtrList<dimensionedScalar> M0(energyGroups);
PtrList<dimensionedScalar> alfaM(energyGroups);
List<scalar> M0_(neutronTransportDict.lookup("M0"));
List<scalar> alfaM_(neutronTransportDict.lookup("alfaM"));

// Transport [1/m]
PtrList<dimensionedScalar> TR0(energyGroups);
PtrList<dimensionedScalar> alfaTR(energyGroups);
List<scalar> TR0_(neutronTransportDict.lookup("TR0"));
List<scalar> alfaTR_(neutronTransportDict.lookup("alfaTR"));

// Sigma_f * nu [1/m]
PtrList<dimensionedScalar> NSF0(energyGroups);
PtrList<dimensionedScalar> alfaNSF(energyGroups);
List<scalar> NSF0_(neutronTransportDict.lookup("NSF0"));
List<scalar> alfaNSF_(neutronTransportDict.lookup("alfaNSF"));

// Sigma_f * Ef [N]
PtrList<dimensionedScalar> SP0(energyGroups);
PtrList<dimensionedScalar> alfaSP(energyGroups);
List<scalar> SP0_(neutronTransportDict.lookup("SP0"));
List<scalar> alfaSP_(neutronTransportDict.lookup("alfaSP"));

// Inverse velocity [s/m]
PtrList<dimensionedScalar> IV0(energyGroups);
List<scalar> IV0_(neutronTransportDict.lookup("IV0"));

// Prompt spectrum [-]
PtrList<dimensionedScalar> chiPro(energyGroups);
List<scalar> chiPro_(neutronTransportDict.lookup("chiPro"));

// Delayed spectrum [-]
PtrList<dimensionedScalar> chiDel(energyGroups);
List<scalar> chiDel_(neutronTransportDict.lookup("chiDel"));

// Scattering (matrix) [1/m]
PtrList<dimensionedScalar> S0(energyGroups*energyGroups);
PtrList<dimensionedScalar> alfaS(energyGroups*energyGroups);
List<scalar> S0_(neutronTransportDict.lookup("S0"));
List<scalar> alfaS_(neutronTransportDict.lookup("alfaS"));

// INITIALIZATION

forAll(M0_, eI)
	{
		M0.set(eI, new dimensionedScalar("M"+Foam::name(eI)+"_0", dimless/dimLength, M0_[eI]));
		TR0.set(eI, new dimensionedScalar("TR"+Foam::name(eI)+"_0", dimless/dimLength, TR0_[eI]));
		NSF0.set(eI, new dimensionedScalar("NSF"+Foam::name(eI)+"_0", dimless/dimLength, NSF0_[eI]));
		SP0.set(eI, new dimensionedScalar("SP"+Foam::name(eI)+"_0", dimensionSet(1,1,-2,0,0,0,0), SP0_[eI]));
		IV0.set(eI, new dimensionedScalar("IV"+Foam::name(eI)+"_0", dimensionSet(0,-1,1,0,0,0,0), IV0_[eI]));
		chiPro.set(eI, new dimensionedScalar("chiPro"+Foam::name(eI), dimless, chiPro_[eI]));
		chiDel.set(eI, new dimensionedScalar("chiDel"+Foam::name(eI), dimless, chiDel_[eI]));

		alfaM.set(eI, new dimensionedScalar("alfa_M"+Foam::name(eI), dimless/dimLength, alfaM_[eI]));
		alfaTR.set(eI, new dimensionedScalar("alfa_TR"+Foam::name(eI), dimless/dimLength, alfaTR_[eI]));
		alfaNSF.set(eI, new dimensionedScalar("alfa_NSF"+Foam::name(eI), dimless/dimLength, alfaNSF_[eI]));
		alfaSP.set(eI, new dimensionedScalar("alfa_SP"+Foam::name(eI), dimensionSet(1,1,-2,0,0,0,0), alfaSP_[eI]));
	}

	forAll(S0_,i)
	{
		S0.set(i, new dimensionedScalar("S"+Foam::name(i)+"_0", dimless/dimLength, S0_[i]));
		alfaS.set(i, new dimensionedScalar("alfa_S"+Foam::name(i), dimless/dimLength, alfaS_[i]));
	}

// Read delayed neutron precursors constants
Info<< "\nReading delayed neutron precursors constants" << endl;

dictionary neutronPrecursorsDict(nuclearProperties.subDict("delayedNeutronPrecursors"));

label nPrecGroups(readInt(neutronPrecursorsDict.lookup("groups")));
PtrList<dimensionedScalar> precLambda(nPrecGroups);
PtrList<dimensionedScalar> precBeta(nPrecGroups);
dimensionedScalar precBetaTot("precBetaTot", dimless, 0);

List<scalar> precLambda_(neutronPrecursorsDict.lookup("lambda"));
List<scalar> precBeta_(neutronPrecursorsDict.lookup("beta"));

if (precLambda_.size() == nPrecGroups && precBeta_.size() == nPrecGroups)
{
    Info<< "    Found " << nPrecGroups << " groups of delayed neutron precursors" << endl;
}
else
{
    FatalError
        << "in " << neutronPrecursorsDict.name() << " , invalid number of delayed neutron precursors constants" << exit(FatalError);
}

forAll(precLambda_, i)
{
    precLambda.set(i, new dimensionedScalar("precLambda" + Foam::name(i + 1), dimless/dimTime, precLambda_[i]));
    precBeta.set(i, new dimensionedScalar("precBeta" + Foam::name(i + 1), dimless, precBeta_[i]));
    precBetaTot += precBeta[i];
}

// Read decay heat precursors constants
Info<< "\nReading decay heat precursors constants" << endl;

dictionary decayHeatPrecursorsDict(nuclearProperties.subDict("decayHeatPrecursors"));

label nDecGroups(readInt(decayHeatPrecursorsDict.lookup("groups")));
PtrList<dimensionedScalar> decLambda(nDecGroups);
PtrList<dimensionedScalar> decBeta(nDecGroups);
dimensionedScalar decBetaTot("decBetaTot", dimless, 0);

List<scalar> decLambda_(decayHeatPrecursorsDict.lookup("lambda"));
List<scalar> decBeta_(decayHeatPrecursorsDict.lookup("beta"));

if (decLambda_.size() == nDecGroups && decBeta_.size() == nDecGroups)
{
    Info<< "    Found " << nDecGroups << " groups of decay heat precursors" << endl;
}
else
{
    FatalError
        << "in " << decayHeatPrecursorsDict.name() << " , invalid number of decay heat precursors constants" << exit(FatalError);
}

forAll(decLambda_, i)
{
    decLambda.set(i, new dimensionedScalar("decLambda" + Foam::name(i + 1), dimless/dimTime, decLambda_[i]));
    decBeta.set(i, new dimensionedScalar("decBeta" + Foam::name(i + 1), dimless, decBeta_[i]));
    decBetaTot += decBeta[i];
}
