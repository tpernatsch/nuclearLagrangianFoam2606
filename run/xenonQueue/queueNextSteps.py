# Next xenon runs (runQueue.py), in order; run all in serial: ./runQueue.sh 8 queueNextSteps.py
# Wall times on 8 local cores (~1.4 s/step at 91k cells, ~3 s/step at 180-216k; 6 cores measured).
#
#   static_*        verification V2 with the final settings (done)    ~1 h in total
#   kin_fb_k1       does 1 smoothing pass keep free fall with the fix  ~0.5 h
#   drop_thesis_*   falling drop to 15 ms, thesis mesh (180k cells)    2.5-6 h each
#   drop_cubic_*    the same on the cubic 15 um mesh (832k cells)      ~30 h
#
# Capillary time-step limit sqrt(rho_mean dx^3 / (2 pi sigma)): 5.8 / 3.2 / 2.1 us for dx = 30 / 20 / 15 um,
# 2.2 us on the thesis mesh (smallest spacing 15.9 um).
# dropForceBalance (per drop zero net CSF - grad p_rgh): without it the drop is held at ~0.01 g t (kin_base),
# with it 0.99 g t (kin_fb). xenonDroplet has it on by default.

def static(res, dt):
    # static drop, Dliq x1000 (verification), cut-cell start, one curvature-smoothing pass, capillary time step
    return dict(name="static_smooth_k1_" + res, case="xenonDropStatic", args=["-res", res, "-init", "smooth"],
                set={"system/controlDict": {"maxDeltaT": dt}, "system/fvSolution": {"nCurvatureSmooth": 1}},
                note="V2 with final settings")


def final(nCS=0, **files):
    # cut-cell start, capillary time step, drop force balance, momentum budget log
    s = {"system/controlDict": {"maxDeltaT": "2e-6"},
         "system/fvSolution": {"nCurvatureSmooth": nCS, "dropForceBalance": "liquid",
                               "dropMomentumBudget": "liquid"}}
    s.update(files)
    return s


RUNS = [
    static("coarse", "5e-6"),
    static("medium", "3e-6"),
    static("fine", "2e-6"),

    # kin_fb (queueKinematics.py) + 1 smoothing pass (without the fix it made the fall 1.5x too fast, kin_nf_k1)
    dict(name="kin_fb_k1", case="xenonDropStatic", args=["-res", "medium", "-init", "smooth"],
         set={"constant/g": {"value": "(0 0 -9.81)"},
              "system/controlDict": {"endTime": "4e-3", "maxDeltaT": "3e-6", "writeInterval": "5e-4"},
              "system/fvSolution": {"nCurvatureSmooth": 1, "dropForceBalance": "liquid",
                                    "dropMomentumBudget": "liquid"}},
         note="drop force balance + 1 curvature-smoothing pass, 4 ms fall"),

    # falling drop, nominal Dliq 2.006e-9, Y0 4e-8, Newman closure
    dict(name="drop_thesis_baseline", case="xenonDroplet", args=[],
         set={"system/fvSolution": {"dropForceBalance": "none"}},
         note="thesis settings: staircase start, stock curvature, maxDeltaT 5 us, no force balance"),
    dict(name="drop_thesis_fb", case="xenonDroplet", args=["-init", "smooth"], set=final(),
         note="cut-cell start, capillary time step, drop force balance"),
    dict(name="drop_thesis_fb_k1", case="xenonDroplet", args=["-init", "smooth"], set=final(1),
         note="drop_thesis_fb + 1 curvature-smoothing pass"),
    dict(name="drop_thesis_fb_Dl5e-10", case="xenonDroplet", args=["-init", "smooth"],
         set=final(**{"constant/phaseProperties": {"Dliq": "5e-10"}}),
         note="drop_thesis_fb, Stokes-Einstein Dliq (removal ~ sqrt(Dliq))"),

    dict(name="drop_cubic_fb", case="xenonDroplet", args=["-cubic", "-init", "smooth"], set=final(),
         note="mesh check: drop_thesis_fb on the cubic 15 um core"),
]

# Gravity in the force balance (dropForceGravity): the reconstructed gravity term alone gave 0.75 g on the
# graded thesis mesh (drop_thesis_fb, 0.88 g t at 15 ms) and 1.07 g on the cubic one (drop_cubic_fb)
G = {"dropForceGravity": "true"}
RUNS += [
    dict(name="kin_fbg", case="xenonDropStatic", args=["-res", "medium", "-init", "smooth"],
         set={"constant/g": {"value": "(0 0 -9.81)"},
              "system/controlDict": {"endTime": "4e-3", "maxDeltaT": "3e-6", "writeInterval": "5e-4"},
              "system/fvSolution": {"dropForceBalance": "liquid", "dropMomentumBudget": "liquid", **G}},
         note="kin_fb with gravity in the force balance, 4 ms fall"),
    dict(name="drop_thesis_fbg", case="xenonDroplet", args=["-init", "smooth"],
         set=final(**{"system/fvSolution": {**final()["system/fvSolution"], **G}}),
         note="drop_thesis_fb with gravity in the force balance"),
    dict(name="drop_cubic_fbg", case="xenonDroplet", args=["-cubic", "-init", "smooth"],
         set=final(**{"system/fvSolution": {**final()["system/fvSolution"], **G}}),
         note="drop_cubic_fb with gravity in the force balance"),
]
