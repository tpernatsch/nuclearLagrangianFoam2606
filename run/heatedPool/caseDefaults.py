"""
caseDefaults.py - suggested answers for setupCase.py, plus the settings that are not asked.

Edit this file to change what the questionnaire suggests, or run
    python3 setupCase.py --defaults
to build a case from it directly.  SI units throughout (m, s, K, W).
"""

# ----------------------------------------------------------------------------- case
GEOMETRY = "A"                 # "A" annulus (5 / 8.6 mm x 75 mm) | "B" cylinder (7 mm x 35 mm)
CASE_NAME = "{geo}_{mesh}{bl}_{mode}"   # placeholders: {geo} {mesh} {bl} {mode}
OUTPUT_DIR = "cases"           # relative to this file

# ----------------------------------------------------------------------------- geometry
A_INNER_DIAMETER = 5.0e-3
A_OUTER_DIAMETER = 8.6e-3
A_HEIGHT = 75.0e-3

B_DIAMETER = 7.0e-3
B_HEIGHT = 35.0e-3
B_CORE_HALF_SIDE = 0.50        # O-grid core half-side / R
B_CORE_ARC_RADIUS = 0.567      # core-edge arc midpoint radius / R  (= B_CORE_HALF_SIDE -> straight edges)

# ----------------------------------------------------------------------------- mesh
MESH_LEVEL = "medium"          # coarse | medium | fine | xfine | xxfine  (cell sizes in setupCase.MESH_LEVELS)
BOUNDARY_LAYER = True          # True: geometric wall refinement on no-slip / heat-transfer walls
BL_FIRST_CELL = 30e-6          # first-cell thickness at refined walls [m]
BL_GROWTH = 1.15               # max cell-to-cell expansion ratio inside the refined layer
MESH_TYPE = "3D"               # "3D" | "wedge" (axisymmetric, one cell around: the steady flows are axisymmetric)
WEDGE_ANGLE = 5.0              # [deg] wedge opening (wedge patches front/back; B: axis as a collapsed edge)

# ----------------------------------------------------------------------------- liquid (salt)
# User data:  rho = 5.12 - 7.98e-4*T [g/cm3]   mu = 0.307*exp(21138/(R*T)) [mPa s]
#             cp = 924 J/kg/K                  k = 1.16 W/m/K
# rho, cp, k, beta are taken at T_REF = 946.15 K (Boussinesq: re-evaluate if T_REF changes).
RHO = 4364.9                   # 5.12 - 7.98e-4*946.15 = 4.3650 g/cm3
CP = 924.0
K = 1.16
BETA = 1.828e-4                # 7.98e-4 / 4.3650  [1/K]
# Viscosity: "arrhenius" -> mu(T) = mu0*exp(Ea/(R T)) (OpenFOAM ArrheniusNewtonian, tangent at T_REF,
# within ~1 % over T_REF +-60 K; k is kept constant by a coded fvOption) | "constant" -> MU
VISCOSITY = dict(model="arrhenius", mu0=0.307e-3, Ea=21138.0)    # mu0 [Pa s], Ea [J/mol]
MU = 4.51e-3                   # used only with model "constant": mu(946.15 K) = 4.51 mPa s

# ----------------------------------------------------------------------------- heat source
SALT_POWER = 20.0              # heat generation per salt mass [W/g] -> x rho = 87.3 W/cm3 at T_REF
POWER_FOLLOWS_DENSITY = True   # q_vol = SALT_POWER x rho(T): T source = SALT_POWER/cp (1 - beta (T - T_REF))
                               # (False: SALT_POWER x rho(T_REF) everywhere)
Q_VOL = 1.0e8                  # [W/m3] used only when SALT_POWER = None

# ----------------------------------------------------------------------------- temperatures
T_INIT = 946.15                # initial liquid temperature [K]  (673 C = mean of measured 710/690/620 C)
T_REF = 946.15                 # Boussinesq reference temperature [K]  (~ expected mean)

# ----------------------------------------------------------------------------- cladding (Hastelloy N)
# k(T) from the material table (linear interpolation, held outside 200-700 C); emissivity 0.25 +- 0.03.
CLADDING = dict(
    k=22.6,                    # [W/m/K] single value (thin shell, lumped estimates): ~670 C
    k_table=[(473.15, 13.1), (573.15, 14.4), (673.15, 16.5), (773.15, 18.0), (873.15, 20.3), (973.15, 23.6)],
    rho=8860.0,
    cp=578.0,                  # [J/kg/K] PLACEHOLDER (literature value at ~700 C; affects transients only)
    generation=2.0,            # [W/g] heat generated in the cladding (walls, inner ring, bottom plate)
    emissivity=0.25,
    outer_diameter=17.0e-3,    # A and B: 17 mm cladding OD
    bottom_thickness=5.0e-3,   # A and B: 5 mm bottom plate
    inner_ring_diameter=2.5e-3,  # A: inner cladding ring, ID 2.5 mm (OD = liquid ID 5 mm)
)
# Not modelled yet: the cladding above the liquid level.

# Thermocouples [C], steady state, three plateaus (the readings step up together over time).
#   PIN 1 = A: T-1/T-2/T-3 on the axis (inside the inner ring) at the top/middle/bottom;
#              T-4..6, T-7..9, T-10..12 on the cladding outer surface at the bottom/middle/top, different theta
#   PIN 2 = B: T-18 outer surface bottom, T-19/T-20 middle (different theta), T-23 top
TC_READINGS = [
    {"T-1": 711, "T-2": 690, "T-3": 612, "T-4": 554, "T-5": 560, "T-6": 548, "T-7": 572, "T-8": 551, "T-9": 571,
     "T-10": 471, "T-11": 540, "T-12": 570, "T-18": 632, "T-19": 617, "T-20": 618, "T-23": 603},
    {"T-1": 721, "T-2": 700, "T-3": 621, "T-4": 564, "T-5": 570, "T-6": 557, "T-7": 582, "T-8": 562, "T-9": 581,
     "T-10": 482, "T-11": 550, "T-12": 580, "T-18": 642, "T-19": 627, "T-20": 629, "T-23": 613},
    {"T-1": 728, "T-2": 706, "T-3": 628, "T-4": 571, "T-5": 577, "T-6": 565, "T-7": 589, "T-8": 569, "T-9": 589,
     "T-10": 491, "T-11": 559, "T-12": 588, "T-18": 649, "T-19": 635, "T-20": 636, "T-23": 621},
]
TC_PLATEAU = 0                 # which plateau (0, 1, 2) matches the simulated power
TC_EXCLUDE = ["T-10"]          # 70-100 K below T-11/T-12 at the same height: cold spot or detached
TC_LAYOUT = {                  # (height fraction, [thermocouples averaged over theta])
    "A": dict(axis=[(0.0, ["T-3"]), (0.5, ["T-2"]), (1.0, ["T-1"])],
              outer=[(0.0, ["T-4", "T-5", "T-6"]), (0.5, ["T-7", "T-8", "T-9"]), (1.0, ["T-10", "T-11", "T-12"])]),
    "B": dict(outer=[(0.0, ["T-18"]), (0.5, ["T-19", "T-20"]), (1.0, ["T-23"])]),
}


def _profiles():
    heights = {"A": 75.0e-3, "B": 35.0e-3}
    out = {}
    for geo, kinds in TC_LAYOUT.items():
        out[geo] = {}
        for kind, pts in kinds.items():
            zs, Ts = [], []
            for f, names in pts:
                vals = [TC_READINGS[TC_PLATEAU][n] for n in names if n not in TC_EXCLUDE]
                zs.append(f * heights[geo])
                Ts.append(sum(vals) / len(vals) + 273.15)
            out[geo][kind] = dict(z=zs, T=Ts)
    return out


THERMOCOUPLES = _profiles()    # {"A": {"axis": {z, T}, "outer": {z, T}}, "B": {"outer": {z, T}}}  [m], [K]
TC_PROFILE = "clamped_top"

# Wall model:
#   "1D"    : conduction through the wall thickness only (axial cladding conduction < ~3 %, see
#             claddingAnalysis.py). Needed by the modes "cladding" and "ring".
#   "shell" : thin-shell finite-area model (regionFaModels thermalShell) on the cladding walls, modes
#             "coefficient" / "flux" / "adiabatic" only (no cladding heat generation yet).
WALL_MODEL = "1D"

# ----------------------------------------------------------------------------- thermal boundaries
# mode "split"       : the heat leaving the salt (generation + ring heat) is shared by mean-flux ratios
#                      (HEAT_SPLIT) between outerWall (reference shape along z), bottom and freeSurface
#                      (uniform). Flux-only: the level comes from LEVEL_PIN.
# mode "cladding"    : coded Robin BC q = H (T_w - T_ref(z)) from steady radial conduction with heat
#                      generation through the Hastelloy cladding (outer wall: tube to the 17 mm OD,
#                      bottom: 5 mm plate), k(T) from CLADDING k_table. Outside of the cladding:
#                        outside = "measured": T_o(z) = thermocouple profile on the outer surface
#                        outside = "gas"     : convection h + radiation (emissivity) to Ta
#                      h = "auto" (gas): calibrated so the energy balance gives T_target at
#                        target = "innerWall" (A axis thermocouples, mean) | "liquid" | "wall" | "surface"
#                      dT_internal = T at the target point - area-mean T of the losing walls
#                      (Tw_loss in postProcessing/heatBalance_*.dat): 0 first, then from a run.
#                      h = "<patch>": same gas (h, Ta) as that patch.
# mode "ring"        : A inner ring (closed inside): its own heat flows into the salt as a fixed flux;
#                      impose_tc = True imposes the axis thermocouple profile instead (fixed T).
# mode "coefficient" : externalWallHeatFluxTemperature: layers + gas h + radiation, no generation
# mode "relative"    : heat flux out = ratio x the current mean heat flux out of patch "of" (per area)
# mode "adiabatic" | "flux" (q_loss [W/m2], > 0 = out) | "balance" (1D only, constant flux, no feedback)
# Any entry can be given per geometry as {"A": ..., "B": ...}. innerWall only exists for geometry A.
THERMAL_BC = {
    "outerWall":   dict(mode="cladding", outside="measured"),   # outer-surface thermocouples (theta mean)
    "innerWall":   dict(mode="ring", impose_tc=True),           # A: fixed at the axis thermocouples (- ring offset)
    "bottom":      dict(mode="relative", of="outerWall", ratio=0.34),  # mean flux 0.34 x side (reference, per area)
    "freeSurface": dict(mode="relative", of="outerWall", ratio=0.17),  # 0.5 x bottom = 0.17 x side (user data)
}
# Alternatives: heat split (all modes "split", level from LEVEL_PIN) or cladding + gas with h calibrated:
#   "outerWall": dict(mode="cladding", outside="gas", Ta=573.15, emissivity=CLADDING["emissivity"], h="auto",
#                     target={"A": "innerWall", "B": "surface"}, T_target={"A": 946.48, "B": 900.0},
#                     dT_internal={"A": 0.0, "B": 0.0}),
#   "bottom":    dict(mode="cladding", outside="gas", h="outerWall", emissivity=CLADDING["emissivity"]),
#   "freeSurface": dict(mode="adiabatic"),

# Heat split (mode "split"): mean flux on the bottom / on the side wall = 0.34 (SALIENT-03 reference model,
# referenceProfiles.py), top free surface / bottom = 0.5 (user data); side shape from the reference model.
HEAT_SPLIT = dict(bottom_over_side=0.34, top_over_bottom=0.5, shape="analysis/reference_shapes.csv",
                  ref_height=36e-3,      # reference capsule height [m]
                  top_zone=3.6e-3)       # top zone of the reference shape (its top 10 %), kept at this absolute
                                         # thickness on our wall (cooling by the cladding above the liquid)

# Level of T with flux-only boundaries: "innerWallTC" (A: inner-wall mean = axis thermocouple spline mean
# - ring offset) | "volume" (volume mean = T_INIT)
LEVEL_PIN = {"A": "innerWallTC", "B": "volume"}
# Wall probes at the thermocouple heights (fractions of the height): A inner wall, B outer wall
TC_HEIGHT_FRACTIONS = {"A": [0.0, 0.5, 1.0], "B": [0.0, 0.5, 1.0]}

# ----------------------------------------------------------------------------- radiation in the liquid
# Participating, grey, non-scattering liquid; walls opaque, grey, diffuse.
#   model "none" | "P1" (cheap, needs a*L >> 1) | "fvDOM" (any optical thickness, ~10x cost)
#   absorption : grey absorption coefficient a [1/m] of the liquid  (PLACEHOLDER: the key input)
#   emissivity : of each surface as seen from the liquid
#   nPhi/nTheta: fvDOM rays per octant (4*nPhi*nTheta rays in 3D); solverFreq: solve every N iterations
# OpenFOAM emits sigma*T^4 (refractive index 1): a liquid with index n emits n^2 times more.
RADIATION = dict(
    model="none",
    absorption=100.0,
    emissivity={"outerWall": 0.25, "innerWall": 0.25, "bottom": 0.25, "freeSurface": 1.0},  # Hastelloy N walls
    nPhi=2, nTheta=2,
)

# ----------------------------------------------------------------------------- run control
RUN_MODE = "steadyThenTransient"   # steady | transient | steadyThenTransient
# Start from a finished case (path relative to this file, or absolute): the internal fields of its latest
# time go into 0.orig, the boundary conditions stay this case's. Same mesh required (same blockMeshDict);
# a decomposed case must be reconstructed first. None -> uniform T_INIT, at rest.
INITIAL_FIELDS = None
STEADY_ITERATIONS = 4000
STEADY_WRITE_INTERVAL = 500
TRANSIENT_END_TIME = 60.0      # [s]  (after the steady stage; for a cold start use several thermal time constants)
TRANSIENT_WRITE_INTERVAL = 5.0 # [s]
MAX_CO = 0.8
MAX_DELTA_T = 0.05             # [s]
HEAT_BALANCE_EVERY = 10        # write the heat balance every N time steps / iterations
N_PROCS = 4

# ----------------------------------------------------------------------------- SLURM cluster
# enabled: write Allrun.slurm (submit from the case directory: sbatch Allrun.slurm) instead of Allrun.
#   ntasks : #SBATCH -n (replaces N_PROCS; numberOfSubdomains is set from SLURM_NTASKS at run time).
#            Wedge cases stay serial (1 task).
#   account: #SBATCH -A ("" -> no line);  extra: more #SBATCH options, e.g. ["-t 48:00:00", "-p cpu"]
CLUSTER = dict(enabled=False, ntasks=64, account="eu", extra=[])

# ----------------------------------------------------------------------------- Lagrangian fission-product particles
# One-way coupled NuclearCloud (src/myIntermediate) riding on the transient stage through the fvOption
# nuclearCloudSourceSingle (src/nuclearCloudFvOptions): no momentum/heat back to the salt, no dispersion,
# no packing, no decay heat, no isotropy, no damping.
#   diameters      : one manualInjection per diameter (injectorID 0, 1, ...), all parcels released at SOI
#   n_per_injection: parcels per injector (nParticle 1: one particle per parcel)
#   release_duration: > 0 -> each injector releases its positions one by one (file order, i.e. random
#                    in space) evenly over [SOI, SOI + release_duration] (manualInjection 'duration')
#   SOI            : release time [s] after the start of the transient stage (the cloud is only built
#                    at the 2nd time step, so it must be > 0)
#   margin         : free fraction of the gap (A) / radius (B) and of the height kept next to every surface;
#                    positions random (seeded), uniform in the 3D volume (r drawn with pdf ~ r), all unique.
#                    Wedge: angle within 90 % of the wedge; the parcels move across it (Brownian, lift)
#                    and the wedge faces reflect them.
#   injection      : "domain"   -> random unique points in the whole liquid (0.25 mm off the walls), manualInjection
#                    "central"  -> random unique points in the central region (margin, below)
#                    "cellZone" -> cellZoneInjection over the cellZone 'liquid' (whole mesh, made by topoSet in
#                                  Allrun): random points in every cell, ~n_per_injection, all at SOI (no duration)
#   walls          : interaction on every wall and the free surface; the wedge faces rebound elastically,
#                    which is exactly the axisymmetric equivalent of 100 % recycling through them
#   particle       : constantProperties (mu0 = salt mu(TRef) unless given)
PARTICLES = dict(
    enabled=False,
    diameters=[1e-9, 1e-8, 1e-7, 1e-6, 1e-5],   # [m]
    n_per_injection=10000,
    SOI=1.0,                                    # [s]
    release_duration=2.0,                       # [s] positions released one by one over this time (0 = all at SOI)
    injection="domain",     # "domain": whole liquid, "central": region with margin, "cellZone": cellZoneInjection
    margin=0.25,            # central region only
    seed=1,
    walls="stick",
    brownian="BrownianMotion",          # stock Li-Ahmadi (exact displacement for dt >> tauP) or "BoundedBrownianMotion" (see README)
    brownian_lambda=1e-12,
    u_integration="analytical",         # exact for frozen forces at any dt/tauP; "CrankNicolson" does not damp for dt >> tauP (README)
    lift="McLaughlinMei",               # shear lift: "McLaughlinMei" (Saffman x McLaughlin J(eps)/J_inf, Mei's fit; small Re),
                                        # "SaffmanMei" (Saffman x Mei finite-Re factor; needs eps >> 1) or "none"
    max_co=0.5,
    log_frequency=1000,
    force_magnitudes=True,      # cloud function ForceMagnitudes (needs it registered for nuclear parcels)
    particle=dict(rho0=10980.0, T0=883.0, Cp0=296.0, kappa0=115.0, mu0=None),
)

# ----------------------------------------------------------------------------- probes
# Points at angle 0 (x axis): every radius x every height. Radii outside the liquid are dropped.
# None -> automatic: mid-gap (A) or axis + mid-radius + near wall (B).
PROBE_RADII = None
PROBE_HEIGHT_FRACTIONS = [0.02, 0.1, 0.25, 0.5, 0.75, 0.9, 0.98]
