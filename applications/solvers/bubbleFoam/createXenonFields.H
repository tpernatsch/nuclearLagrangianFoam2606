//Reading Xenon dictionary

Info << "\nReading Xenon properties" << endl;

dictionary xenonDict(nuclearProperties.subDict("xenon"));
dimensionedScalar lambda_xe(xenonDict.lookup("lambda_xe"));
dimensionedScalar MW_xe(xenonDict.lookup("MW_xe"));
dimensionedScalar y_xe(xenonDict.lookup("y_xe"));
dimensionedScalar N_Av(xenonDict.lookup("N_Av"));

// Xenon capture [m^2]
PtrList<dimensionedScalar> alfaXE(energyGroups);
List<scalar> alfaXE_(neutronTransportDict.lookup("alfaXE"));

forAll(alfaXE_,eI)
{
  alfaXE.set(eI, new dimensionedScalar("alfa_XE"+Foam::name(eI), dimensionSet(0,2,0,0,0,0,0), alfaXE_[eI]));
}

//adding fissionRate for Xe production
          volScalarField fissionRate
          (
              IOobject
              (
                  "fissionRate",
                  runTime.timeName(),
                  mesh,
                  IOobject::READ_IF_PRESENT,
                  IOobject::AUTO_WRITE
              ),
              mesh
          );

//adding Xenon mass, for each phase and total one
          volScalarField xenonSolute
          (
              IOobject
              (
                  "xenonSolute",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ,
                  IOobject::AUTO_WRITE
              ),
              rho2*phase2.Y("xenon")*alpha2
          );

          volScalarField xenonGas
          (
              IOobject
              (
                  "xenonGas",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ,
                  IOobject::AUTO_WRITE
              ),
              rho1*phase1.Y("xenon")*alpha1
          );

          volScalarField xenonMass
          (
              IOobject
              (
                  "xenonMass",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ,
                  IOobject::AUTO_WRITE
              ),
              xenonSolute + xenonGas
          );

//adding Xenon absorption cross sections
PtrList<volScalarField> XE(energyGroups);

forAll(XE,eI)
{
    XE.set
    (
        eI,
        new volScalarField
        (
            IOobject
            (
                "XE"+Foam::name(eI),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            xenonMass*(N_Av/MW_xe)*alfaXE[0]
        )
    );
}

for (int eI=0; eI<energyGroups; ++eI)
{
  XE[eI] = xenonMass*(N_Av/MW_xe) * alfaXE[eI];
}

//adding another field for Xe capture
volScalarField xeCapture
(
    IOobject
    (
        "xeCapture",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    0*alfaXE[0]*flux[0]
);
