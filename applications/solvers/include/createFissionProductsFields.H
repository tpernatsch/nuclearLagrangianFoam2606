Info<< "\nReading fission products fields (" << nFPSpecies << ")" << endl;

// Set fission products fields size
Dfp.setSize(nFPSpecies);
fp.setSize(nFPSpecies);
fps.setSize(nFPSpecies);
fpp.setSize(nFPSpecies);
fpd.setSize(nFPSpecies);
d.setSize(nFPSpecies);
qdecay_wall.setSize(nFPSpecies);
gamma.setSize(nFPSpecies);

forAll(fp, i)
{
    // Diffusion coefficient
    Dfp.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Dfp" + name(i + 1),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Dfp" + name(i + 1), dimArea/dimTime, 0)
        )
    );
    
    Info<< " \n Reading field fp" << i + 1 << endl;
    
    // Fission products concentration
    fp.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "fp" + name(i + 1),
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    Info<< " \n Reading field fps" << i + 1 << endl;

    //Dissolved fission products concentration
    fps.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "fp" + name(i + 1) + "s",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("", dimless/dimVolume, 0)
        )
    );

    Info<< " \n Reading field fpp" << i + 1 << endl;

    // Precipited fission products concentration
    fpp.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "fp" + name(i + 1) + "p",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("", dimless/dimVolume, 0)
        )
    );

    Info<< " \n Reading field fpd" << i + 1 << endl;

    // Fission products deposition rate
    fpd.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "fp" + name(i + 1) + "d",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("", dimless/dimArea/dimTime, 0)
        )
    );

    // Deposition velocity (output only)
    gamma.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "gamma" + name(i + 1),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("gamma" + name(i + 1), dimVelocity, 0)
        )
    );


    // Total fission products deposited
    d.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "d" + name(i + 1),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("d" + name(i + 1), dimless/dimArea, 0)
        )
    );
    
    // Decay power of deposited fission products
    qdecay_wall.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "qdecay_wall" + name(i + 1),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("qdecay_wall" + name(i + 1), dimPower/dimArea, 0)
        )
    );
}
