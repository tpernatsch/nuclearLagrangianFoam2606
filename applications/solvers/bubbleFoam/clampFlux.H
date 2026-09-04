// Numerical safety net for bubbleFoam's transient mode.
//
// Unlike the reference designs this solver's neutronics were adapted from
// (msfrSimpleFoam's updateKeff.H/normFlux.H, and msfrParcelFoam's
// criticality-search path), bubbleFoam's transient mode (criticalityCalculation=0)
// holds Keff fixed with NO active per-step correction -- intentionally, so
// that a real reactivity transient (e.g. from gas injection) is visible
// rather than being renormalized away every timestep like a steady-state
// SIMPLE solve would. That means nothing here stops the flux solve from
// numerically diverging if a timestep is locally too large for the flux
// equation's own (very stiff, IV-weighted) dynamics.
//
// This block only engages as a guard rail: if prompt power has diverged
// far outside any physically sane bound, rescale flux (and everything
// linear in it) back down before the corrupted power reaches the
// temperature equation. Under normal operation -- including a real,
// legitimate power excursion from gas injection -- this never fires.

{
    // Check LOCAL maxima, not just the domain-integrated total -- the
    // failure mode seen throughout this case's debugging has consistently
    // been a single cell (or a handful) reaching an extreme value while
    // the domain average looks only mildly elevated, so a
    // domain-integrated check alone would miss exactly the blowup this
    // guard is meant to catch.
    scalar promptPower = fvc::domainIntegrate(promptHeatSource).value();

    // Absolute local-flux bound, sized off this case's own initial uniform
    // flux guess (5e14, see regen_neutronics_fields.py / 0/flux0..5) with a
    // generous margin. Not a general-purpose/portable constant -- if this
    // case's flux scale changes substantially, this bound should move with
    // it.
    scalar fluxBound = 5e14*50.0;

    scalar maxLocalFlux = 0;
    forAll(flux, eI)
    {
        maxLocalFlux = max(maxLocalFlux, max(flux[eI]).value());
    }

    scalar boundHi = 10.0*nominalPower;

    bool blownUp =
        !std::isfinite(promptPower)
     || mag(promptPower) > boundHi
     || !std::isfinite(maxLocalFlux)
     || maxLocalFlux > fluxBound;

    if (blownUp)
    {
        WarningInFunction
            << "Neutron flux diverged (prompt power = " << promptPower
            << " W, nominal = " << nominalPower << " W; max local flux = "
            << maxLocalFlux << ", bound = " << fluxBound
            << ") -- clamping flux back within bounds."
            << " This is a numerical safety net, not real physics; if it"
            << " fires repeatedly, the timestep needs to be smaller, not"
            << " this clamp." << endl;

        scalar normFactor =
        (
            std::isfinite(maxLocalFlux) && maxLocalFlux > fluxBound
          ? maxLocalFlux/fluxBound
          : (
                std::isfinite(promptPower) && mag(promptPower) > VSMALL
              ? promptPower/nominalPower
              : 1e6
            )
        );

        forAll(flux, eI)
        {
            flux[eI] = flux[eI]/normFactor;
        }

        NSFflu = NSFflu/normFactor;
        SPflu = SPflu/normFactor;
        promptHeatSource = promptHeatSource/normFactor;
        fissionRate = fissionRate/normFactor;
        xeCapture = xeCapture/normFactor;

        forAll(flux, eI)
        {
            scaNeu[eI] = scaNeu[eI]/normFactor;
            othGroFisNeu[eI] = othGroFisNeu[eI]/normFactor;
        }
    }
}
