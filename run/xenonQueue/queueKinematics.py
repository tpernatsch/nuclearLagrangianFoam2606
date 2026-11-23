# Drop kinematics (runQueue.py): the falling drop should be in free fall (drag time ~0.4 s), U = g t.
# Medium static mesh (20 um cubic) with gravity, 4 ms (g t = 0.039 m/s, fall 0.08 mm), cut-cell start,
# capillary time step. ~0.5 h each on 8 cores.

def fall(name, note, fvSolution=None, args=()):
    fv = {"nCurvatureSmooth": 0}
    fv.update(fvSolution or {})
    return dict(name="kin_" + name, case="xenonDropStatic",
                args=["-res", "medium", "-init", "smooth", *args],
                set={"constant/g": {"value": "(0 0 -9.81)"},
                     "system/controlDict": {"endTime": "4e-3", "maxDeltaT": "3e-6", "writeInterval": "5e-4"},
                     "system/fvSolution": fv},
                note=note)


RUNS = [
    fall("base", "3 alpha sub-cycles, stock curvature"),
    fall("nsub1", "1 alpha sub-cycle", {"nAlphaSubCycles": 1}),
]

# Fix: zero each drop's discrete net CSF force (solver PIMPLE csfZeroNetForce)
NF = {"csfZeroNetForce": "liquid"}
RUNS += [
    fall("nf", "zero net CSF force", NF),
    dict(name="static_nf_medium", case="xenonDropStatic", args=["-res", "medium", "-init", "smooth"],
         set={"system/controlDict": {"maxDeltaT": "3e-6"},
              "system/fvSolution": {"nCurvatureSmooth": 1, **NF}},
         note="static_smooth_k1_medium + zero net CSF force (spurious currents)"),
]
RUNS += [
    # free fall without surface tension: is the late overshoot of kin_nf CSF-related?
    fall("sigma0", "no surface tension, 4 ms", None) | {"set": {**fall("x", "")["set"],
                                                              "constant/phaseProperties": {"sigma": 0}}},
]
RUNS += [fall("nf_k1", "zero net CSF force + 1 curvature-smoothing pass", {**NF, "nCurvatureSmooth": 1})]

# Momentum budget per drop (dropMomentumBudget log) and the mass-consistency Sp term
BUD = {**NF, "dropMomentumBudget": "liquid"}
RUNS += [
    fall("nf_bud", "zero net CSF force, momentum budget log", BUD),
    fall("nf_sp", "zero net CSF force + Sp term, momentum budget log", {**BUD, "momentumSp": "true"}),
]

# Combined fix: per drop, net (CSF - grad p_rgh) = gas-side pressure force, in the U reconstruction
FB = {"dropForceBalance": "liquid", "dropMomentumBudget": "liquid"}
RUNS += [
    fall("fb", "drop force balance (CSF + pressure), budget log", FB),
    dict(name="static_fb_medium", case="xenonDropStatic", args=["-res", "medium", "-init", "smooth"],
         set={"system/controlDict": {"maxDeltaT": "3e-6"},
              "system/fvSolution": {"nCurvatureSmooth": 1, **FB}},
         note="static_smooth_k1_medium + drop force balance (spurious currents)"),
    fall("fb_nf", "drop force balance + zero net CSF force, budget log", {**FB, **NF}),
]
