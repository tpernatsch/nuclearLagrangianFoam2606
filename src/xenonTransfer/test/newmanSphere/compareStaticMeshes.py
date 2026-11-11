#!/usr/bin/env python3
"""V2: static drop on several uniform meshes (run/xenonDropStatic).

Usage: compareStaticMeshes.py <Dliq> <label>=<caseDir> [<label>=<caseDir> ...]
Per mesh: Newman error of the liquid inventory (checkNewmanCase), drop age,
total xenon, initial diameter, and the parasitic velocity from the log
(uniform mesh: Umax = Co_max*dx/dt).
"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import checkNewmanCase as cn


def courant(case, dx):
    # Co is printed before setDeltaT, i.e. with the previous step's dt
    out, co, dtPrev = [], None, None
    for line in open(os.path.join(case, "log.xenonTransferInterFoam")):
        if line.startswith("Courant Number mean"):
            co = float(line.split("max:")[1])
        elif line.startswith("deltaT = "):
            dt = float(line.split()[2])
            if co is not None and dtPrev:
                out.append(co*dx/dtPrev)
            co, dtPrev = None, dt
    return out


def dxOf(case):
    s = open(os.path.join(case, "system", "blockMeshDict")).read()
    return 1e-3*float(re.search(r"^dx\s+([0-9.eE+-]+);", s, re.M).group(1))


def main():
    Dliq = float(sys.argv[1])
    rows = []
    for arg in sys.argv[2:]:
        label, case = arg.split("=", 1)
        ts, ages, (M0, d0) = cn.parseLog(os.path.join(case, "log.xenonTransferInterFoam"))
        liq = cn.series(case, "ZXeMass_in_liquid")
        tot = cn.series(case, "ZXeMass_total")
        dx = dxOf(case)
        err, last = [], None
        for t in sorted(liq):
            i = min(range(len(ts)), key=lambda k: abs(ts[k] - t))
            Fo = 4*Dliq*ages[i]/d0**2
            R = cn.remaining(Fo)
            err.append(liq[t]/M0/R - 1)
            last = (t, ages[i]/t, Fo, 1 - liq[t]/M0, 1 - R, tot[t]/M0)
        U = courant(case, dx)
        n = len(U)
        # Boundary xenon (boundaryXe: XEOUT t rate integral; outflow > 0)
        out = 0.0
        for line in open(os.path.join(case, "log.xenonTransferInterFoam")):
            if line.startswith("XEOUT "):
                out = float(line.split()[3])
        tEnd = sorted(tot)[-1]
        imbalance = (tot[tEnd] + out - M0)/M0
        rows.append((label, dx, d0, max(err, key=abs), last, max(U),
                     sum(U[n//2:])/max(len(U[n//2:]), 1), out/M0, imbalance))

    print("%-7s %-6s %-8s %-10s %-8s %-8s %-9s %-10s %-10s %-10s %-7s %-7s" %
          ("mesh", "dx um", "d0 um", "max err", "theta/t", "F CFD", "F Newman", "Mtot/M0",
           "out/M0", "imbalance", "Umax", "U late"))
    for (label, dx, d0, e, (t, a, Fo, F, FN, mt), umax, ulate, out, imb) in rows:
        print("%-7s %-6.1f %-8.2f %+-10.2e %-8.5f %-8.5f %-9.5f %-10.7f %+-10.2e %+-10.2e %-7.3g %-7.3g" %
              (label, dx*1e6, d0*1e6, e, a, F, FN, mt, out, imb, umax, ulate))


if __name__ == "__main__":
    main()
