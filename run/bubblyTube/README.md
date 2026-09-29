# bubblyTube — Euler-Euler bubbly flow in an open 30 mm tube (OpenFOAM v2606)

Case generator for a vertical tube (D = 30 mm) with a free surface and a centred gas injector.
It is set up after F.L. de Wit's MSc thesis (TU Delft, 2025), referred to here as **[dW]**, with the baffle removed.

| File | Purpose |
|---|---|
|  | **Run this.** It asks the questions (numbered options, suggested answer in [brackets], Enter accepts it), then writes one complete case per flowrate into  |
|  | The suggested answers, plus the settings that are not asked (e.g. the liquid level, fixed at 300 mm). Normally not edited |
|  | Reduces a finished run to CSV tables for comparison with LDA and imaging |



The questionnaire has six steps (plus optional advanced settings):

1. **Operating conditions:** liquid (AG40 / water); flowrate(s) in L/min, format , suggested range 0.1-0.5, comma-separated for a sweep; flowmeter reference conditions
2. **Geometry and mesh:** tube height, format  (the liquid level stays at 0.30 m); mesh level, with cell counts and cell/d<sub>b</sub> shown
3. **Bubbles:** diameter from the thesis table or a typed value ( mm); size model
4. **Physics:** turbulence (laminar / LaheyKEpsilon / kOmegaSSTSato / mixtureKEpsilon), drag, lift and its wall damping, wall lubrication, turbulent dispersion and C<sub>td</sub> (skipped when laminar)
5. **Measurements:** LDA plane heights ( m) and the LDA point layout
6. **Run control:** simulated time, start of averaging, processors

Advanced (optional): time step, Courant number, relaxation, correctors, faceMomentum, aspect ratio, C<sub>vm</sub>, write interval, probe frequency, force output, orifice diameter.

Wrong inputs are rejected with a message. Values outside the suggested range ask for confirmation.
Every case stores its answers in :



## What gets generated

- **Mesh:** a curved O-grid, as in [dW] App. C.3. It has a square core with arc edges (`CORE_ARC_RADIUS`) and 4 outer blocks, all hexahedral. z points up, with the origin at the centre of the bottom. The patches are `bottom` and `wall` (walls) and `atmosphere` (open top).
  The `medium` level reproduces the quality metrics of [dW] Table C.1 exactly (max aspect ratio 2.76, non-orthogonality 22.6°, skewness 0.52).
- **Injection:** mass and energy sources (`scalarSemiImplicitSource`) on a small cellZone above the orifice, as in [dW] fvOptions. The mass flow is ṁ = ρ<sub>ref</sub>Q, with the reference conditions set by `FLOW_REFERENCE`.
- **Initial state and outlet:** the liquid column is set to `LIQUID_HEIGHT` with `setFields`, with air above it. The open top is a `prghPressure` + `inletOutlet` free surface, equivalent to [dW] BC III.
- **Run-time outputs, for each height in `LDA_HEIGHTS`:**
  - probes at every LDA point (radii × angles);
  - two mean-profile diameters (`ldaLines`) and a full plane (`ldaPlanes`);
  - `voidZi`: the instantaneous area-averaged void fraction;
  - `meanVoidZi`: the area average of α<sub>mean</sub>;
  - `flowZi`: ∫U<sub>z,mean</sub><0 dA, ∫U<sub>z,mean</sub>>0 dA and A<sub>down</sub>. This is the open-tube analogue of the [dW] downcomer flow rate.

  Plus the holdup, the gas outflow, vertical α<sub>water</sub> lines for the free-surface level, and the time averages (`fieldAverage`) from `AVERAGE_START` onwards.
- **`Allrun`:** runs blockMesh → checkMesh → topoSet → setFields → decomposePar → solver → reconstructPar (latest time).

## Where this deviates from [dW], and why

| Item | [dW] (OpenFOAM-8, Foundation) | Here (v2606, ESI) |
|---|---|---|
| Solver | `multiphaseEulerFoam` | `reactingTwoPhaseEulerFoam`. The v2606 `multiphaseEulerFoam` is a different code with only drag and heat-transfer closures (no lift, turbulent dispersion or BIT turbulence) |
| Liquid thermo | `eConst` + `rPolynomial` | `hConst` + `perfectFluid` (constant density at 1 atm). [dW] gives no rPolynomial coefficients for AG40 |
| Lift wall damping | `wallDamped`, `zeroWallDist = 0.0002` | `wallDamped` + `cosine`, `Cd = 1`. This is the same function as [dW] Eq. 8.13 with d<sub>range</sub> = d<sub>b</sub>; `zeroWallDist` does not exist in v2606 |
| Turbulent dispersion | LopezDeBertodano Ctd 2.0 (setup B), Burns Ctd 1, σ 0.7 (setup C) | Default LopezDeBertodano Ctd 2.0. **Burns diverged in tests** (AG40, d<sub>b</sub> 6.75 mm), both from rest and when switched on at t = 0.2 s. The option is kept |
| Time step | fixed 1 ms, maxCo 5 | adaptive, maxCo 0.5, maxΔt 1 ms (≈1 ms in practice). `FIXED_TIME_STEP = 1e-3` reproduces [dW] |
| Relaxation | 0.85 | 1.0 (`RELAXATION`) |
| Bubble size | constant 3 mm (downcomer-representative) | Riser sizes of [dW] §5.2, 5 → 8.5 mm for 0.1 → 0.5 L/min. Without a baffle there is no downcomer, so the rising bubbles drive the flow |
| Domain | 0.30 m (BC I/II) or 0.50 m (BC III, "to keep it stable") | `TUBE_HEIGHT` = 0.35 m. If the free surface becomes unstable, increase it |

**The v2606 turbulent-dispersion models scale with the liquid k or ν<sub>t</sub>.** With `TURBULENCE = "laminar"` there is therefore no turbulent dispersion at all.

## Bubble size vs mesh — the numbers to discuss

The bubbles are 5–8.5 mm across in a 30 mm tube (d<sub>b</sub>/D = 0.17–0.28). No mesh can have cells both larger than the bubble (the Euler-Euler averaging assumption) and small enough to resolve the recirculation.

| Level | N<sub>c</sub>/N<sub>r</sub>/Δz | cells (0.35 m) | Δ core / arc / Δz [mm] | h [mm] | Δ/d<sub>b</sub> at 5 mm | at 8.5 mm |
|---|---|---|---|---|---|---|
| coarse | 8/4/3.0 mm | 22.5 k | 1.9 / 2.9 / 3.0 | 2.2 | 0.37–0.6 | 0.22–0.35 |
| medium | 12/4/2.0 mm | 58.8 k | 1.25 / 2.0 / 2.0 | 1.6 | 0.25–0.4 | 0.15–0.24 |
| fine | 16/6/1.5 mm | 149 k | 0.94 / 1.5 / 1.5 | 1.2 | 0.19–0.3 | 0.11–0.18 |
| extrafine | 22/8/1.15 mm | 360 k | 0.68 / 1.07 / 1.15 | 0.9 | 0.14–0.23 | 0.08–0.14 |

Points relevant to the choice:
1. **Mesh convergence in [dW] depended on the turbulence model.** With AG40, Lahey k-ε was *not* mesh-converged (GCI 70–127 %, y<sup>+</sup> in the buffer layer). k-ω SST Sato was converged at the medium mesh (GCI<sub>21</sub> = 0.54 %).
2. **The Tomiyama lift changes sign at d<sub>b</sub> = 6.26 mm** for AG40, which falls inside the measured range. Lift pushes bubbles to the wall at 0.1–0.2 L/min and towards the centre at 0.3 L/min and above. The bubble size therefore matters more than the mesh.
3. **Cosine wall damping acts over Cd·d<sub>b</sub> = 5–8.5 mm**, so lift is reduced over 1/3 to more than 1/2 of the radius.

## Things to check

- `BUBBLE_TABLE` (in `caseDefaults.py`): [dW] only states the end points (5 mm and 8.5 mm). The intermediate values are interpolated linearly, so replace them with values read off [dW] Fig. 5.2 or with your own imaging.
- `TUBE_HEIGHT`: 0.35 m (your spec) vs 0.45 m ([dW] rig).
- `FLOW_REFERENCE`: [dW] converts L/min with ρ = 1.204 kg/m³ (20 °C). A Bronkhorst MFC is usually calibrated in *normal* L/min (0 °C).
- Probe files grow with `len(LDA_HEIGHTS) × len(LDA_ANGLES) × len(LDA_RADII)` points written every `PROBE_WRITE_EVERY` steps after `AVERAGE_START`. The default (400 points, every step, 30 s) produces a few hundred MB per field. Choose 4 radial paths, or a larger probe interval in the advanced settings, if that is too much.
