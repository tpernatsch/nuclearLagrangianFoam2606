# Next xenon runs (runQueue.py), cheapest first. Wall times: measured cost on 6 local cores
# (~1.4 s/step at 91k cells, ~3 s/step at 180-216k); the 32/64-core figures assume ~70 % parallel efficiency.
#
#   static_*        verification V2 with the final settings          ~1 h local in total, ~20 min on 32 cores
#   drop_thesis_*   falling drop to 15 ms, thesis mesh (180k cells)  2.5-6 h local each,  ~45 min-2 h on 32 cores
#   drop_cubic_*    the same on the cubic 15 um mesh (832k cells)    ~30 h local,          ~4-5 h on 64 cores
#
# Capillary time-step limit sqrt(rho_mean dx^3 / (2 pi sigma)): 5.8 / 3.2 / 2.1 us for dx = 30 / 20 / 15 um,
# 2.2 us on the thesis mesh (smallest spacing 15.9 um).

def static(res, dt):
    # static drop, Dliq x1000 (verification), cut-cell start, one curvature-smoothing pass, capillary time step
    return dict(name="static_smooth_k1_" + res, case="xenonDropStatic", args=["-res", res, "-init", "smooth"],
                set={"system/controlDict": {"maxDeltaT": dt}, "system/fvSolution": {"nCurvatureSmooth": 1}},
                note="V2 with final settings")


FINAL = {"system/controlDict": {"maxDeltaT": "2e-6"}, "system/fvSolution": {"nCurvatureSmooth": 1}}

RUNS = [
    static("coarse", "5e-6"),
    static("medium", "3e-6"),
    static("fine", "2e-6"),

    # falling drop, nominal Dliq 2.006e-9, Y0 4e-8, Newman closure
    dict(name="drop_thesis_baseline", case="xenonDroplet", args=[],
         note="current settings: staircase start, stock curvature, maxDeltaT 5 us"),
    dict(name="drop_thesis_final", case="xenonDroplet", args=["-init", "smooth"], set=FINAL,
         note="cut-cell start, 1 smoothing pass, capillary time step"),
    dict(name="drop_thesis_final_Dl5e-10", case="xenonDroplet", args=["-init", "smooth"],
         set=dict(FINAL, **{"constant/phaseProperties": {"Dliq": "5e-10"}}),
         note="Stokes-Einstein Dliq (removal ~ sqrt(Dliq))"),

    dict(name="drop_cubic_final", case="xenonDroplet", args=["-cubic", "-init", "smooth"], set=FINAL,
         note="mesh check: cubic 15 um core"),
]
