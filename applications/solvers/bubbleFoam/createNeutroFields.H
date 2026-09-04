//Reading cross section dictionary
Info<< "\n Reading nuclear properties \n" << endl;

IOdictionary nuclearProperties
(
    IOobject
    (
        "nuclearProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar TrefXS(nuclearProperties.lookup("TrefXS"));
dimensionedScalar rhoReference(nuclearProperties.lookup("rhoReference"));
dimensionedScalar nubar(nuclearProperties.lookup("nubar"));

volScalarField rhok=alpha2*rho2/rhoReference;
volScalarField rhod=(alpha1*rho1 + alpha2*rho2)/rhoReference;
volScalarField Tfuel=thermo2.T();

scalar neutronTransport = runTime.controlDict().lookupOrDefault<scalar>("neutronTransport", 0);
scalar criticalityCalculation = runTime.controlDict().lookupOrDefault<scalar>("criticalityCalculation", 0);
scalar THCoupling = runTime.controlDict().lookupOrDefault<scalar>("THCoupling", 1);

  #include "readNuclearProperties.H"

// Create albedo flux and diffCoef

    Info<< "Reading field fluxAlbedo\n" << endl;
    volScalarField fluxAlbedo
    (
        IOobject
        (
            "fluxAlbedo",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//initialize variables for the albedo routine
    volScalarField nullFlux
    (
        IOobject
        (
            "nullFlux",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field diffCoeffAlbedo\n" << endl;
    volScalarField diffCoeffAlbedo
    (
        IOobject
        (
            "diffCoeffAlbedo",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

	PtrList<volScalarField> flux(energyGroups);
	forAll(flux, eI)
	{
		flux.set
		(
			 eI,
			 new volScalarField
			 (
          IOobject
          (
              "flux"+Foam::name(eI),
              runTime.timeName(),
              mesh,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE
          ),
          mesh
			 )
		);
	}

	PtrList<volScalarField> secondFlux(energyGroups);
	forAll(secondFlux, eI)
	{
		secondFlux.set
		(
			 eI,
			 new volScalarField
			 (
          IOobject
          (
              "secondFlux"+Foam::name(eI),
              runTime.timeName(),
              mesh,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE
          ),
          mesh
			 )
		);
	}

//declare Keff
// Matches msfrParcelFoam's convention (KeffInit read from controlDict)
// -- also needed so a criticality-search result can be carried into a
// transient run instead of every run silently assuming Keff=1.
scalar Keff = runTime.controlDict().lookupOrDefault<scalar>("Keff", 1);
scalar Keff_new=Keff;
scalar res = 1;
scalar reactivity = 0;

scalar voidManual = runTime.controlDict().lookupOrDefault<scalar>("voidFraction", 1);
scalar pumpControl = runTime.controlDict().lookupOrDefault<scalar>("pumpControl", 0);
scalar hexControl = runTime.controlDict().lookupOrDefault<scalar>("hexControl", 1);
scalar reactorVolume = runTime.controlDict().lookupOrDefault<scalar>("reactorVolume", 1);
scalar pumpHeigth = runTime.controlDict().lookupOrDefault<scalar>("pumpHeigth", 1);
scalar externalCycles = runTime.controlDict().lookupOrDefault<scalar>("externalCycles", 1);
scalar KeffControl = runTime.controlDict().lookupOrDefault<scalar>("KeffControl", 1);
//scalar friction = runTime.controlDict().lookupOrDefault<scalar>("friction", 1);

// Delayed neutron precursor and decay heat precursor constants (precLambda,
// precBeta, precBetaTot, decLambda, decBeta, decBetaTot) are read from the
// nuclearProperties subDicts in readNuclearProperties.H

    Info<< "Delayed neutron precursor (" << nPrecGroups << ") prec\n" << endl;

    PtrList<volScalarField> prec(nPrecGroups);

    forAll(prec, i)
    {
        Info<< "    Reading field prec" << i + 1 << endl;
        prec.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "prec" + Foam::name(i + 1),
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }

//Read Schmidt number
    Info<< "\n Reading Schmidt Number for neutron precursors \n" << endl;

    IOdictionary precursorTransportProperties
    (
      IOobject
      (
        "precursorTransportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
      )
    );

	dimensionedScalar Sc(precursorTransportProperties.lookup("Sc"));
	dimensionedScalar Sct(precursorTransportProperties.lookup("Sct"));

//read decay heat precursor

    Info<< "Decay heat precursor (" << nDecGroups << ") dec\n" << endl;

    PtrList<volScalarField> dec(nDecGroups);

    forAll(dec, i)
    {
        Info<< "    Reading field dec" << i + 1 << endl;
        dec.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "dec" + Foam::name(i + 1),
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }

//declare XS fields

	volScalarField logT("logT", log(Tfuel/TrefXS));

    PtrList<volScalarField> NSF(energyGroups);
    PtrList<volScalarField> SP(energyGroups);

    forAll(NSF,eI)
    {
        NSF.set
        (
            eI,
            new volScalarField
            (
                IOobject
                (
                    "NSF"+Foam::name(eI),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                0*(NSF0[0] + alfaNSF[0] * logT)*rhok
            )
        );

        SP.set
        (
            eI,
            new volScalarField
            (
                IOobject
                (
                    "SP"+Foam::name(eI),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                0*(SP0[0] + alfaSP[0] * logT)*rhok
            )
        );
    }

	List<PtrList<volScalarField > > S(energyGroups);
	forAll(S,eI)
	{
		S[eI].resize(energyGroups);

		forAll(S[eI],eJ)
		{
			S[eI].set
			(
				eJ,
				new volScalarField
				(
					IOobject
					(
						"S"+Foam::name(eI)+Foam::name(eJ),
						runTime.timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					0*(S0[0] + alfaS[0] * logT)*rhok
				)
			);
		}
	}

    PtrList<volScalarField> M(energyGroups);
    PtrList<volScalarField> TR(energyGroups);
    PtrList<volScalarField> IV(energyGroups);
    PtrList<volScalarField> F(energyGroups);

    forAll(M,eI)
    {
        M.set
        (
            eI,
            new volScalarField
            (
                IOobject
                (
                    "M"+Foam::name(eI),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                0*(M0[0] + alfaM[0] * logT)*rhok
            )
        );

        TR.set
        (
            eI,
            new volScalarField
            (
                IOobject
                (
                    "TR"+Foam::name(eI),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                0*(TR0[0] + alfaTR[0] * logT)*rhok
            )
        );

        IV.set
        (
            eI,
            new volScalarField
            (
                IOobject
                (
                    "IV"+Foam::name(eI),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                0*IV0[0]
            )
        );

        F.set
        (
            eI,
            new volScalarField
            (
                IOobject
                (
                    "F"+Foam::name(eI),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                0*NSF[0] * (1-precBetaTot) * chiPro[0]
            )
        );
    }

	volScalarField NSFflu("NSFflu", 0*NSF[0] * flux[0]);
	volScalarField SPflu("SPflu", 0*SP[0] * flux[0]);

//delayed neutrons source
// scattered neutrons
    PtrList<volScalarField> scaNeu(energyGroups);
    PtrList<volScalarField> othGroFisNeu(energyGroups);

    forAll(scaNeu,eI)
    {
        scaNeu.set
        (
            eI,
            new volScalarField
            (
                IOobject
                (
                    "scaNeu"+Foam::name(eI),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                0*S[0][0] * flux[0]
            )
        );

        othGroFisNeu.set
        (
            eI,
            new volScalarField
            (
                IOobject
                (
                    "othGroFisNeu"+Foam::name(eI),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                0*(1-precBetaTot) * NSF[0] * flux[0] * chiPro[0]
            )
        );
    }

    for (int eI=0; eI<energyGroups ; ++eI)
    {
    		SP[eI] = (SP0[eI] + alfaSP[eI] * logT)*rhok;
    		NSF[eI] = (NSF0[eI] + alfaNSF[eI] * logT)*rhok;
    		F[eI] = NSF[eI] * (1.0-precBetaTot) * chiPro[eI];
    		M[eI] = (M0[eI] + alfaM[eI] * logT)*rhok;
    		TR[eI] = (TR0[eI] + alfaTR[eI] * logT)*rhok;
        IV[eI] = IV0[eI];

        for (int eJ=0; eJ<energyGroups ; ++eJ)
    		{
    				int eIJ = energyGroups*eI + eJ;
    				S[eI][eJ] = (S0[eIJ] + alfaS[eIJ] * logT)*rhok;
		    }

	      NSFflu = NSFflu + NSF[eI] * flux[eI];
	      SPflu  = SPflu  + SP[eI] * flux[eI];
  }


	volScalarField NSFflu_old("NSFflu_old", NSFflu);
  scalar FS_integral = fvc::domainIntegrate(NSFflu).value();
  scalar FS_integral_old = fvc::domainIntegrate(NSFflu_old).value();

    PtrList<volScalarField> delNeuSource(energyGroups);
    forAll(delNeuSource,eI)
    {
        delNeuSource.set
        (
            eI,
            new volScalarField
            (
                IOobject
                (
                    "delNeuSource"+Foam::name(eI),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                prec[0]*precLambda[0]*0//*rho2*alpha2
            )
        );
    }

	for (int eI=0; eI<energyGroups ; ++eI)
	{
		delNeuSource[eI] = chiDel[eI] * precLambda[0] * prec[0] * 0;
		forAll(prec, j)
		{
			delNeuSource[eI] += chiDel[eI] * precLambda[j] * prec[j];
		}

		for (int eJ=0; eJ<energyGroups ; ++eJ)
		{
			if (eI!=eJ)
  			{
				   scaNeu[eI] = scaNeu[eI] + S[eJ][eI] * flux[eJ];
           othGroFisNeu[eI] = othGroFisNeu[eI] + (1-precBetaTot) * NSF[eJ] * flux[eJ] * chiPro[eI];
        }
    }
	}

    PtrList<volScalarField> delNeuSource_old(energyGroups);
    PtrList<volScalarField> scaNeu_old(energyGroups);
    PtrList<volScalarField> othGroFisNeu_old(energyGroups);

    forAll(delNeuSource_old,eI)
    {
        delNeuSource_old.set
        (
            eI,
            new volScalarField
            (
                IOobject
                (
                    "delNeuSource_old"+Foam::name(eI),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                delNeuSource[0]
            )
        );

        scaNeu_old.set
        (
            eI,
            new volScalarField
            (
                IOobject
                (
                    "scaNeu_old"+Foam::name(eI),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                scaNeu[0]
            )
        );

        othGroFisNeu_old.set
        (
            eI,
            new volScalarField
            (
                IOobject
                (
                    "othGroFisNeu_old"+Foam::name(eI),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                othGroFisNeu[0]
            )
        );
    }

    volScalarField promptHeatSource
    (
        IOobject
        (
            "promptHeatSource",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        ((1-decBetaTot)*SPflu)
    );

    volScalarField decayHeatSource
    (
        IOobject
        (
            "decayHeatSource",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dec[0]*decLambda[0]*0 //*rho2*alpha2 //modify here
    );

    forAll(dec, i)
    {
        decayHeatSource += decLambda[i]*dec[i];
    }

// other declaration
Info<< "Declare all\n" << endl;

//create power density field

    Info<< "Power density field\n" << endl;


    volScalarField powDens
    (
        IOobject
        (
            "powDens",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        promptHeatSource + decayHeatSource
    );

//create total power

    Info<< "creating totalPower\n" << endl;

    scalar parTotalVolume = 0;

    forAll(mesh.cells(),celli)
    {
        parTotalVolume += mesh.V()[celli];
    }

    reduce(parTotalVolume, sumOp<scalar>());

//create diffusion coefficient of precursors

    volScalarField diffEff = (thermo2.nu()/Sc + phase2.nut()/Sct);//*rho2*alpha2;

    Info << nl << nl << "totalPower = " << fvc::domainIntegrate(powDens).value() << " W" << nl << nl << endl;
    Info << nl << nl << "parTotalVolume = " << parTotalVolume << nl << nl << endl;

    scalar nominalPower = runTime.controlDict().lookupOrDefault<scalar>("nominalPower", 1e9);

// Create volPow

    Info<< "volPow \n" << endl;

    Info<< "Reading field volPow\n" << endl;
    volScalarField volPow
    (
        IOobject
        (
            "volPow",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
