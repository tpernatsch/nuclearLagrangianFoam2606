"""
C_inputs.py - CONFIDENTIAL input data for heatedPool geometry "C" (graphite crucible 2 of the
predecessor irradiation). Source: the irradiation PIE report supplied by the user (confidential).

NEVER commit this file or copy its numbers into tracked files, commit messages, slides or web
queries. The folder run/heatedPool/confidential/ is gitignored; generated cases (cases/) are too.
SI units unless stated.  "fig." values are read off plots (+- a few %).
"""

# ----------------------------------------------------------------------------- salt
SALT = dict(
    composition="78LiF-22ThF4",
    T_melt=568.0 + 273.15,           # [K]
    mass=6.1344e-3,                  # [kg]  crucible 2 sample
    th232_mass=3.5564e-3,            # [kg]
)
# Liquid properties, 78LiF-22ThF4: Benes & Konings (2012), Comprehensive Nuclear Materials 3.13, Table 4
# (public literature, not from the report):
#   rho = 5543 - 1.25 T [kg/m3] (Porter & Meaker, eqn 7);  mu = 0.365 exp(2735/T) [mPa s] (Chervinskij, eqn 8)
#   cp = 1.0 J/g/K (estimate, ideal + 11 %);  k ~ 1.5 W/m/K (estimate at 1023 K);  T_melt 841 K
T_REF_C = 995.0                     # [K] salt volume mean of the EOI run with the reference wall (994.3 K)
SALT_PROPS = dict(
    rho_a=5543.0, rho_b=1.25,        # rho = a - b T
    mu0=0.365e-3, Ea=2735.0 * 8.314462618,   # mu = mu0 exp(Ea/(R T)), Ea = 22740 J/mol
    cp=1000.0, k=1.5,
)
#  T [K]   rho     beta [1/K]   mu [mPa s]   H [mm] (mass / rho / bore area)
#   908   4408.0   2.84e-4       7.42        36.16
#   960   4343.0   2.88e-4       6.30        36.70
#   980   4318.0   2.89e-4       ~5.95       ~36.9
#  1000   4293.0   2.91e-4       5.62        37.13
# Room-temperature 36.4 mm -> apparent solid density 4379 kg/m3 (voids / meniscus): use the mass-based height.

# ----------------------------------------------------------------------------- geometry (crucible 2)
COLUMN_DIAMETER = 7.0e-3            # salt bore diameter
COLUMN_HEIGHT_RT = 36.4e-3          # measured at room temperature (meniscus convex before irradiation)
# Operating height = mass / (rho_liquid(T) * pi D^2/4)  (graphite expansion neglected) -> once rho(T) known
FREE_SURFACE = "flat"               # user decision 2026-11-16 (post-irradiation the salt crept up the walls)

CRUCIBLE = dict(
    material="graphite T-950",
    outer_diameter=26.0e-3,          # ~25-26.5 mm (drawing of the same crucible design)
    length=75.0e-3,
    bore_depth_from_top=60.0e-3,     # top 15 mm: wider plenum (~18.4 mm) + lid thread; below: salt bore
    bottom_thickness=15.0e-3,        # 75 - 60
    tc_circle_diameter=21.5e-3,      # thermocouple holes (1.1 mm) pitch circle
    monitor_circle_diameter=23.1e-3, # fluence monitor sets
    k=None,                          # T-950 conductivity: NOT AVAILABLE -> generic irradiated value + sensitivity
    generation=3.0,                  # [W/g] nuclear heating of the graphite (reference CFD assumption)
)

# ----------------------------------------------------------------------------- operating point
# Ru-103 half-life 39.3 d -> the scanned Ru-103 reflects the last ~2-3 months: model end of irradiation.
SALT_POWER_EOI = 26.6               # [W/g] crucible 2, last cycle (26.9 -> 26.6); total ~163 W
POWER_HISTORY_NOTE = "8.2 W/g at start -> 26.6 W/g at EOI; drop ~day 345 (moved to a lower-flux position)"
TC_WALL_EOI = 635.0 + 273.15        # [K] crucible 2 thermocouples (in the graphite wall, salt mid-height),
                                    # histogram peak ~635 C; component table 620-645 C (fig./table)
TC_HEIGHT_FRACTION = 0.5            # only one axial position
# Axial wall profile (user decision 2026-11-16): B's outer thermocouple profile (bottom / middle / top at
# height fractions 0, 0.5, 1; same TC_PLATEAU and "clamped_top" spline), shifted so that the middle
# point equals TC_WALL_EOI; applied at the thermocouple circle (r = tc_circle_diameter / 2).
TC_PROFILE_SOURCE = "B shifted to the mid-height reading"

GRAPHITE = dict(
    k=30.0,          # [W/m/K] PLACEHOLDER: irradiated nuclear graphite at ~900 K (generic); sensitivity 20-60
    rho=1860.0,      # [kg/m3] reference CFD
    cp=710.0,        # [J/kg/K] reference CFD (transients only)
    emissivity=0.8,  # generic graphite (only used with an outside = "gas" boundary)
)


# Wall model of C (user decision 2026-11-16):
#   "reference": the salt-wall T of the reference CFD base case (REF_CFD_BASE wall_inner, its full model of the
#                graphite, gas gaps and containments), shifted by T_tc(C) - T_tc(ref) and imposed on the salt wall
#                (a 'cladding' 0.01 % of R thick: Dirichlet). Valid because crucible 2 at EOI has the base case's
#                heat per length within 2 % (26.6 x 7.0^2 vs 23 x 7.6^2). Validated: with it imposed our salt model
#                gives the reference's U (13.1 vs 12-12.5 mm/s), shear (0.36 vs 0.35 Pa), Tmax - TC (168 vs 165 K).
#                Free surface adiabatic (reference: stagnant gas above, no radiation).
#   "graphite" : radial conduction through the graphite to the TC circle, B's profile shape shifted (first try;
#                graphiteConduction.py shows the wrong axial trend). Kept as a sensitivity.
WALL = "reference"


def preset(b_outer):
    """Settings overriding caseDefaults.py for geometry C (applied by setupCase.apply_geometry_preset).
    b_outer: B's outer thermocouple profile {z, T} (z at fractions 0, 0.5, 1 of B's height)."""
    import math
    p = SALT_PROPS
    T = T_REF_C
    rho = p["rho_a"] - p["rho_b"] * T
    R = COLUMN_DIAMETER / 2
    H = SALT["mass"] / rho / (math.pi * R * R)
    HB = b_outer["z"][-1]
    shift = TC_WALL_EOI - b_outer["T"][1]                     # B's middle point -> C's reading
    tc = dict(z=[z / HB * H for z in b_outer["z"]], T=[t + shift for t in b_outer["T"]])
    g = GRAPHITE
    extra = {}
    if WALL == "reference":
        rb = REF_CFD_BASE
        sh = TC_WALL_EOI - rb["T_tc"]
        tc = dict(z=[1e-3 * z / (1e3 * rb["H"]) * H for z, _ in rb["wall_inner"]],
                  T=[t + 273.15 + sh for _, t in rb["wall_inner"]])
        extra = dict(TC_PROFILE="pchip")          # monotone between the read-off points
    out = dict(
        C_DIAMETER=COLUMN_DIAMETER, C_HEIGHT=H,
        RHO=rho, CP=p["cp"], K=p["k"], BETA=p["rho_b"] / rho,
        VISCOSITY=dict(model="arrhenius", mu0=p["mu0"], Ea=p["Ea"]),
        MU=p["mu0"] * math.exp(p["Ea"] / (8.314462618 * T)),
        T_REF=T, T_INIT=T,
        SALT_POWER=SALT_POWER_EOI, POWER_FOLLOWS_DENSITY=True,
        CLADDING=dict(k=g["k"], k_table=None, rho=g["rho"], cp=g["cp"], generation=CRUCIBLE["generation"],
                      emissivity=g["emissivity"],
                      outer_diameter=CRUCIBLE["tc_circle_diameter"],   # 'measured' T at the thermocouple circle
                      bottom_thickness=CRUCIBLE["bottom_thickness"],
                      inner_ring_diameter=None),
        THERMAL_BC={
            "outerWall": dict(mode="cladding", outside="measured"),
            "innerWall": dict(mode="adiabatic"),
            "bottom": dict(mode="cladding", outside="measured"),      # floor: T_o = profile at z = 0
            "freeSurface": dict(mode="relative", of="outerWall", ratio=0.17),
        },
        LEVEL_PIN="volume",
        **extra,
    )
    if WALL == "reference":                      # imposed salt-wall T: a vanishing 'cladding', no generation
        out["CLADDING"] = dict(out["CLADDING"], generation=0.0, outer_diameter=COLUMN_DIAMETER * (1 + 1e-4),
                               bottom_thickness=1e-7)
        out["THERMAL_BC"]["freeSurface"] = dict(mode="adiabatic")
    return out, tc

# ----------------------------------------------------------------------------- reference CFD (Fluent, laminar, steady)
# Run for crucibles 1 and 4 (not 2): 23 W/g, thermocouples 590 C.
REF_CFD = dict(
    salt_power=23.0, T_tc=590.0 + 273.15,
    U_up_centre=12e-3, U_up_centre_err=1e-3,      # [m/s] hot salt up in the centre, down at the wall ~same speed
    wall_shear=0.4, wall_shear_err=0.08,          # [Pa] max on the crucible inner surface
    dT_salt_50W=(70.0, 80.0),                     # [K] Tmax - Tmin at ~50 W
    dT_salt_165W=(180.0, 190.0),                  # [K] at 160-165 W
    T_salt_min_vs_tc=10.0,                        # [K] min salt T within 10 K of the thermocouples
)

# Reference CFD report details (CONFIDENTIAL, second report supplied 2026-11-16): Fluent 15, laminar, steady,
# 2D = 3D within 0.1 K. Salt k 1.5, cp 1000, rho(T), mu(T); graphite rho 1860, cp 710, k(T); He-Ne gas
# between salt and filter (70 % Ne); containment gas k tuned to hit the TC temperature. TC at z = 18 mm above
# the salt bottom, r = 10.75 mm. Nuclear heating in all non-salt materials 3 W/g (graphite included).
# Salt column 36 mm high. No radiation mentioned (top of the salt: hot, stably stratified layer).
REF_CFD_BASE = dict(                 # base case, detailed model
    D=7.6e-3, H=36e-3, salt_power=23.0, other_power=3.0, T_tc=590.0 + 273.15, z_tc=18e-3,
    U_up=12e-3, U_down=12e-3, tau_max=0.35,
    # Fig. 4 (read off, +- 2 K): [z mm, T C]
    wall_inner=[(0, 582), (2, 588), (5, 595), (10, 603), (15, 608), (20, 611), (25, 613.5), (30, 614.5),
                (33, 614), (35, 612), (36, 609)],
    axis=[(0, 583), (4, 638), (10, 660), (13, 671), (20, 686), (28, 701), (30, 705), (33, 715), (35, 728),
          (36, 742)],
    graphite_outer=[(-10, 558), (0, 568), (10, 579), (20, 587), (24, 588), (30, 585), (36, 579), (40, 572),
                    (50, 556), (55, 553), (65, 558)],      # r = crucible outer surface
)
# Sensitivity map (Figs. 10-14, 7.6 mm, base settings) at 23 W/g: Tmin - TC ~ -12 K, Tmax - TC ~ 165 K,
# vmax ~ 12.5 mm/s (mu x2: unchanged; k_salt 1: ~10 mm/s), tau_max ~ 0.35 Pa (mu x2: ~0.45).
# Simplified model (D 7.0, no containments/plug/filter): ~40 K hotter everywhere; axis 620 -> 782 C top.
REF_CFD_MAP_23 = dict(dTmin=-12.0, dTmax=165.0, vmax=12.5e-3, tau=0.35)

# ----------------------------------------------------------------------------- Ru-103 gamma scans (crucible 2)
# Horizontal scan across the salt column: collimator 0.5 mm wide x 5 mm high, step 0.5 mm (fig., a.u.)
RU103_HORIZONTAL = dict(
    collimator=(0.5e-3, 5.0e-3),
    wall_peaks=(1.5e5, 1.85e5),      # at x ~ -3 and +3 mm (salt-graphite wall)
    interior=(0.65e5, 0.75e5),       # flat interior
    note="wall / interior ~ 2.2-2.6; deposition of a Ru phase on the graphite wall",
)
# Axial scan of the column: collimator 0.5 mm axial x 15 mm wide, step 0.5 mm (fig., kBq per position)
RU103_AXIAL = dict(
    collimator=(0.5e-3, 15.0e-3),
    top_surface_peak=2.85e6,         # at the free surface
    interior_low=5e4,                # just below the surface
    near_bottom=5e5,                 # sloping rise towards the bottom
    bottom_peak=2.0e6,
    note="bottom rise possibly settling after shutdown, before freezing (report's interpretation)",
)
