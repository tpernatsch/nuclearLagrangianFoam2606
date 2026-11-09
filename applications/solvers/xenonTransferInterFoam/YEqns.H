{
    const label nPhases = fluid.phases().size();

    PtrList<PtrList<volScalarField::Internal>> SusAll(nPhases);
    PtrList<PtrList<volScalarField::Internal>> SpsAll(nPhases);

    // =============================================================
    // PASS 1:
    // Build all Su/Sp source terms first.
    // Important: do NOT solve any Yi/Z equation in this pass.
    // =============================================================
    label phasei = 0;

    for (multiphaseInter::phaseModel& phase : fluid.phases())
    {
        PtrList<volScalarField>& Y = phase.Y();

        SusAll.set
        (
            phasei,
            new PtrList<volScalarField::Internal>(Y.size())
        );

        SpsAll.set
        (
            phasei,
            new PtrList<volScalarField::Internal>(Y.size())
        );

        PtrList<volScalarField::Internal>& Sus = SusAll[phasei];
        PtrList<volScalarField::Internal>& Sps = SpsAll[phasei];

        forAll(Y, i)
        {
            Sus.set
            (
                i,
                new volScalarField::Internal
                (
                    IOobject
                    (
                        "Su" + Y[i].name(),
                        mesh.time().timeName(),
                        mesh
                    ),
                    mesh,
                    dimensionedScalar(dimless/dimTime, Zero)
                )
            );

            Sps.set
            (
                i,
                new volScalarField::Internal
                (
                    IOobject
                    (
                        "Sp" + Y[i].name(),
                        mesh.time().timeName(),
                        mesh
                    ),
                    mesh,
                    dimensionedScalar(dimless/dimTime, Zero)
                )
            );
        }

        forAll(Y, i)
        {
            // Assemble source terms using the same pre-solve fields
            // for all phases.
            fluid.massSpeciesTransfer(phase, Sus[i], Sps[i], Y[i].name());
        }

        ++phasei;
    }

    // =============================================================
    // PASS 2:
    // Now solve Yi/Z after every phase source has already been frozen.
    // =============================================================
    phasei = 0;

    for (multiphaseInter::phaseModel& phase : fluid.phases())
    {
        PtrList<volScalarField>& Y = phase.Y();

        if (!Y.empty())
        {
            phase.solveYi(SusAll[phasei], SpsAll[phasei]);
        }

        ++phasei;
    }
}
