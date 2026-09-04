// Redefinition (already declared)
    rhoCpRef = rhoInfValue*CpRef;

    #include "readNuclearProperties.H"

// Read neutron flux fields
    Info<< "\nReading neutron flux fields (" << nEnergyGroups << ")" << endl;

    PtrList<volScalarField> flux(nEnergyGroups);

    volScalarField fluxTot
    (
        IOobject
        (
            "fluxTot",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless/dimArea/dimTime, 0)
    );

    forAll(flux, i)
    {
        Info<< "    Reading field flux" << i + 1 << endl;
        flux.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "flux" + name(i + 1),
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );

        fluxTot += flux[i]; // update fluxTot with initial flux fields
    }

// Read delayed neutron precursors fields
    Info<< "\nReading delayed neutron precursors fields (" << nPrecGroups << ")" << endl;

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
                    "prec" + name(i + 1),
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }

// Read decay heat precursors fields
    Info<< "\nReading decay heat precursors fields (" << nDecGroups << ")" << endl;

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
                    "dec" + name(i + 1),
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }

// Read fission products fields (if needed)
    PtrList<volScalarField> Dfp;  // Diffusion coefficient
    PtrList<volScalarField> fp;   // Fission products concentration [part./m3]
    PtrList<volScalarField> fps;  // Dissolved fission products concentration
    PtrList<volScalarField> fpp;  // Precipited fission products concentration
    PtrList<volScalarField> fpd;  // Fission products deposition rate [part./m2s]
    PtrList<volScalarField> gamma;  // Fission products deposition velocioty [m/s]
    PtrList<volScalarField> d;    // Total fission products deposited
    PtrList<volScalarField> qdecay_wall;  // Decay heat flux of deposited fission products

    if (fpTransport)
    {
        #include "createFissionProductsFields.H"
    }

// Create neutron transport fields
    Info<< "\nCreating neutron transport constants fields" << endl;

    // Temperature-dependent properties are declared as fields, which populate
    // corresponding pointer lists
    // Dummy initialisation to reference values
    PtrList<volScalarField> D(nEnergyGroups);
    PtrList<volScalarField> Sa(nEnergyGroups);
    PtrList<volScalarField> Sf(nEnergyGroups);
    PtrList<volScalarField> Ss(nEnergyGroups*nEnergyGroups);

    // Net production operators, to be used later
    PtrList<volScalarField> M(nEnergyGroups);

    forAll(D, i)
    {
        word iName = name(i + 1);

        D.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "D" + iName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                D_ref[i]
            )
        );

        Sa.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Sa" + iName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                Sa_ref[i]
            )
        );

        Sf.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Sf" + iName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                Sf_ref[i]
            )
        );

        forAll(D, j)
        {
            label ijInd = nEnergyGroups*i + j;
            word ijName = name(i + 1) + name(j + 1);

            Ss.set
            (
                ijInd,
                new volScalarField
                (
                    IOobject
                    (
                        "Ss" + ijName,
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    Ss_ref[ijInd]
                )
            );
        }

        M.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "M" + iName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                0*Nu[i]*Sf_ref[i]
            )
        );
    }

    volScalarField diffCoeffAlbedo
    (
        IOobject
        (
            "diffCoeffAlbedo",
            runTime.timeName(),
            mesh
        ),
        0*D[0]
    );

// Create other fields

    const volScalarField unitField // unit dimless field to help declaration of other fields
    (
        IOobject
        (
            "unitField",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 1)
    );

    const volScalarField zeroField = 0*unitField; // zero dimless field

    volScalarField Rf // fission rate
    (
        IOobject
        (
            "Rf",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless/dimVolume/dimTime, 0)
    );

    volScalarField q // heat source density
    (
        IOobject
        (
            "q",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimPower/dimVolume, 0)
    );

    volScalarField qPrompt // prompt heat source density
    (
        IOobject
        (
            "qPrompt",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimPower/dimVolume, 0)
    );

    volScalarField qDecay // decay heat source density
    (
        IOobject
        (
            "qDecay",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimPower/dimVolume, 0)
    );

    volScalarField Q // integrated heat source (uniform field)
    (
        IOobject
        (
            "Q",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimPower, 0)
    );

    volScalarField Qprompt // integrated prompt heat source (uniform field)
    (
        IOobject
        (
            "Qprompt",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimPower, 0)
    );

    volScalarField Qdecay // integrated decay heat source (uniform field)
    (
        IOobject
        (
            "Qdecay",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimPower, 0)
    );

 //   volScalarField qPrompt = q;
 //   volScalarField qDecay = q;

    volScalarField NuRf = zeroField*dimensionedScalar("", dimless/dimVolume/dimTime, 0);
    volScalarField EfRf = zeroField*dimensionedScalar("", dimEnergy/dimVolume/dimTime, 0);

    PtrList<volScalarField> neutronSource(nEnergyGroups);
    forAll(neutronSource, i)
    {
        neutronSource.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "neutronSource" + name(i + 1),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("", dimless/dimVolume/dimTime, 0)
            )
        );
    }

    volScalarField rhoNorm
    (
        IOobject
        (
            "rhoNorm",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    );

    volScalarField logT = zeroField;

// Read values related to reactor control
    Info<< "\nReading reactor control parameters" << endl;

    Switch tempCoupling(runTime.controlDict().lookup("tempCoupling"));
    Switch criticality(runTime.controlDict().lookupOrDefault<word>("criticality", "yes"));
    dimensionedScalar KeffInit("KeffInit", dimless, runTime.controlDict().lookupOrDefault<scalar>("Keff", 1));
    dimensionedScalar Qnominal("nominalPower", dimPower, runTime.controlDict());

    volScalarField Keff
    (
        IOobject
        (
            "Keff",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        KeffInit
    );

    volScalarField reactivity
    (
        IOobject
        (
            "reactivity",
            runTime.timeName(),
            mesh
        ),
        (Keff - 1)/Keff
    );

    Info<< "    Nominal reactor power  : " << Qnominal.value()/1E+06 << " MW" << endl;
    Info<< "    Keff                   : " << max(Keff).value() << endl;
    Info<< "    Reactivity             : " << max(reactivity).value()*1e5 << " pcm" << endl;

// Initialise fields
    #include "updateCrossSections.H"
    #include "updateFissionRate.H"
    #include "updatePowerSource.H"
    #include "updateNeutronSource.H"
