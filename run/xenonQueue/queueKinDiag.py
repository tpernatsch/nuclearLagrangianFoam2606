# Drop-kinematics diagnosis (runQueue.py): 0.5 ms toggles of kin_base (queueKinematics.py), which holds the
# drop at 0.07 g t while stock interFoam on the same mesh and fluids reaches 0.84 g t. ~4 min each on 8 cores.

def diag(name, note, extra):
    s = {"constant/g": {"value": "(0 0 -9.81)"},
         "system/controlDict": {"endTime": "5e-4", "maxDeltaT": "3e-6", "writeInterval": "2.5e-4"},
         "system/fvSolution": {"nCurvatureSmooth": 0}}
    for f, kv in extra.items():
        s.setdefault(f, {}).update(kv)
    return dict(name="diag_" + name, case="xenonDropStatic", args=["-res", "medium", "-init", "smooth"],
                set=s, note=note)


RUNS = [
    diag("sigma0", "no surface tension", {"constant/phaseProperties": {"sigma": 0}}),
    diag("noMT", "no mass transfer (C 0)", {"constant/phaseProperties": {"C": 0}}),
    diag("rhog100", "gas 100x denser", {"constant/thermophysicalProperties.gas": {"rho": 0.2615}}),
]

# Bond-number check: the CSF support should fail when the needed curvature asymmetry (Bo/4) grows
RUNS += [
    diag("g10", "10x gravity (Bo 0.14)", {"constant/g": {"value": "(0 0 -98.1)"}}),
    diag("sigma01", "sigma / 10 (Bo 0.14)", {"constant/phaseProperties": {"sigma": 0.02617}}),
]
