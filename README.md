# nuclearLagrangianFoam2606

Personal OpenFOAM user library and solver collection for Lagrangian particle
tracking in nuclear (molten salt reactor) and general multiphase CFD
contexts, built against **OpenFOAM v2606** (OpenCFD/ESI).

## Layout

### Libraries (`src/`)

| Directory | What it is |
| --- | --- |
| `myIntermediate` | Fork of `lagrangian/intermediate`: adds `NuclearCloud`/`NuclearParcel` (decay heat, radiation) and `NuclearEulerCloud`/`NuclearEulerParcel` (two-phase, Euler-Euler-aware particle tracking with a stochastic phase-switch model) on top of the stock `KinematicCloud`/`ThermoCloud` template stack. |
| `myTurbulence` | Fork of the Lagrangian turbulence dispersion/force submodels used by `myIntermediate`'s clouds (RAS dispersion models, two-phase Brownian motion force, ...). |
| `nuclear` | Small support library for nuclear-specific solvers (`msfrParcelFoam`). |
| `nuclearCloudFvOptions` | `fvOption`s that run a `myIntermediate` cloud one-way-coupled on top of an existing solver: `nuclearCloudSource` (single carrier phase) and `nuclearEulerCloudSource` (two named Euler-Euler phases). |
| `kinematicCloudFvOptions` | Single-phase incompressible counterpart of `nuclearCloudFvOptions`, for plain kinematic clouds on top of `pimpleFoam`/`simpleFoam`-like solvers. |
| `flotation` | Mineral-flotation particle-population transport (`flotationSource` fvOption + collision/attachment efficiency, particle velocity/diffusivity and size-distribution submodels), ported from OpenFOAM 13 to v2606. |
| `LogMoM` | Log-normal method-of-moments `diameterModel` for `reactingEuler`/`multiphaseSystem`: tracks a phase's particle size distribution via its scaled number-concentration/surface-area moments instead of a full sectional population balance. |

### Solvers (`applications/solvers/`)

| Solver | What it is |
| --- | --- |
| `nuclearParcelFoam` | Transient buoyant, turbulent, incompressible flow with Lagrangian particle tracking (thermophoresis, decay heat, thermal effects) - typical use: fission-product/noble-metal transport. |
| `msfrParcelFoam` | Same as above, expanded with neutronics and nuclear source terms, fully coupled to flow and temperature - for Molten Salt Reactor (MSFR) applications. |
| `bubbleFoam` | Two compressible fluid phases sharing a common pressure, run-time-selectable phase models. |
| `thermoParcelFoam` | Boussinesq-approximation incompressible PIMPLE solver with a kinematic Lagrangian cloud, surface film modelling, radiation and dynamic mesh support. |

### Utilities (`applications/utilities/`)

| Utility | What it is |
| --- | --- |
| `setLogNormal` | Initialises a phase's `lambda`/`kappai` moment fields from a log-normal size distribution (mean diameter, standard deviation), for use with `LogMoM`. |

## Requirements

- OpenFOAM **v2606** (OpenCFD/ESI distribution), sourced in the shell
  (`source /usr/lib/openfoam/openfoam2606/etc/bashrc` or your install's
  equivalent) so `$WM_PROJECT_DIR`/`$WM_PROJECT_USER_DIR` and the `wmake`
  toolchain are on `PATH`.
- You can clone this repository anywhere - it does not need to sit at
  OpenFOAM's conventional `$HOME/OpenFOAM/$USER-<version>` user directory
  or match the local username.

## Building

From the repository root, with the OpenFOAM environment already sourced,
**every new shell**:

```sh
source etc/sourceme.sh
```

This points `WM_PROJECT_USER_DIR`, `FOAM_USER_APPBIN`/`FOAM_USER_LIBBIN`
and `PATH`/`LD_LIBRARY_PATH` at this repository, wherever it lives -
without it, `wmake`/`wclean` can't find `myIntermediate`/`myTurbulence`
(`-llagrangianMyIntermediate`/`-llagrangianMyTurbulence` "cannot be found"
is this step missing), and built solvers won't be on `PATH`. It must be
*sourced* (`source ...` / `. ...`), not executed - a script can only export
variables into its own subshell, not back into the shell that ran it.

Then build everything:

```sh
./Allwmake
```

This builds every library and solver/utility above, in dependency order
(`myIntermediate` before `myTurbulence` before the two `*CloudFvOptions`
libraries, etc.) into this repository's own `platforms/<arch>/lib` and
`platforms/<arch>/bin`.

To build or clean a single directory by hand instead of the whole repo
(e.g. after editing just one solver), `source etc/sourceme.sh` first, then
`wmake`/`wclean` that directory as usual.

To remove everything `Allwmake` produces (compiled objects, `lnInclude`
link farms and the `platforms/` output directory):

```sh
./Allwclean
```

## Running

Once built, the solvers behave like any other OpenFOAM application: run
them from inside a case directory that provides the usual `system/`,
`constant/` and `0/` setup (e.g. `nuclearParcelFoam` in a case configured
like the stock `reactingParcelFoam` tutorials, plus a `<cloudName>Properties`
dictionary for the Lagrangian cloud). The `*CloudFvOptions` libraries are
instead loaded into an *existing* solver (e.g. `reactingMultiphaseEulerFoam`)
via `libs (...)` in `system/controlDict`, and configured through an entry in
`constant/fvOptions` - see the `Description` block at the top of
`nuclearCloudSource.H`/`nuclearEulerCloudSource.H` for example dictionaries.

No test cases are included in this repository; case setup is expected to
live alongside your own run directories.
