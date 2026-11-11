#!/usr/bin/env python3
"""V-N2: compare the liquid xenon inventory of a gasDiffusionSimple run with
Newman's stagnant sphere, M_l/M0 = 1 - F(Fo), Fo = 4 Dliq theta/d0^2.

Usage: checkNewmanCase.py <caseDir> <Dliq> [out.csv]
Needs the case's ZXeMass_in_liquid / ZXeMass_total volIntegrate functions and
debugDrops (M0 = first MxeDrop, d0 = first Dvol, theta(t) = drop age).
MxeDrop itself is not used: it counts only drop-cluster cells, and liquid
smeared out of the cluster leaves it without being transferred.
"""
import bisect
import glob
import math
import re
import sys

PI = math.pi


def remaining(Fo):
    if Fo <= 0:
        return 1.0
    if Fo < 0.05:
        return 1 - 6*math.sqrt(Fo/PI) + 3*Fo
    return 6/PI**2*sum(math.exp(-n*n*PI*PI*Fo)/(n*n) for n in range(1, 60))


def parseLog(path):
    ts, ages, first = [], [], None
    t = None
    with open(path) as f:
        for line in f:
            if line.startswith("Time = "):
                t = float(line.split()[2])
            elif line.startswith("  drop 0:") and t is not None:
                g = lambda k: float(re.search(r" %s=([-+0-9.eE]+)" % k, line).group(1))
                if first is None:
                    first = (g("MxeDrop"), g("Dvol"))
                if ts and ts[-1] == t:
                    ages[-1] = g("age")
                else:
                    ts.append(t)
                    ages.append(g("age"))
    return ts, ages, first


def series(case, name):
    rows = {}
    for p in sorted(glob.glob("%s/postProcessing/%s/*/volFieldValue*.dat" % (case, name))):
        for line in open(p):
            if line.startswith("#") or not line.strip():
                continue
            a = line.split()
            rows[float(a[0])] = float(a[1])
    return rows


def main():
    case, Dliq = sys.argv[1], float(sys.argv[2])
    out = sys.argv[3] if len(sys.argv) > 3 else None

    ts, ages, (M0, d0) = parseLog("%s/log.xenonTransferInterFoam" % case)
    liq = series(case, "ZXeMass_in_liquid")
    tot = series(case, "ZXeMass_total")

    def age(t):
        i = min(max(bisect.bisect_left(ts, t), 1), len(ts) - 1)
        t0, t1 = ts[i - 1], ts[i]
        return ages[i - 1] + (ages[i] - ages[i - 1])*(t - t0)/max(t1 - t0, 1e-300)

    res = []
    for t in sorted(liq):
        th = age(t)
        Fo = 4*Dliq*th/d0**2
        R = remaining(Fo)
        res.append((t, th, Fo, liq[t]/M0, R, tot.get(t, float("nan"))/M0))

    print("M0=%.6e kg  d0=%.6e m  Dliq=%.3e" % (M0, d0, Dliq))
    print("%-10s %-10s %-10s %-10s %-10s %-10s %-10s" %
          ("t", "theta/t", "Fo", "Ml/M0", "1-F", "relErr", "Mtot/M0"))
    step = max(1, len(res)//20)
    for i, (t, th, Fo, m, R, mt) in enumerate(res):
        if i % step == 0 or i == len(res) - 1:
            print("%-10.4e %-10.6f %-10.4e %-10.6f %-10.6f %-+10.2e %-10.6f" %
                  (t, th/t, Fo, m, R, m/R - 1, mt))
    err = [abs(m/R - 1) for (_, _, _, m, R, _) in res]
    print("samples=%d  max|Ml/M0/(1-F) - 1|=%.3e  final=%.3e  final Mtot/M0=%.6f" %
          (len(res), max(err), err[-1], res[-1][5]))

    if out:
        with open(out, "w") as f:
            f.write("t,theta,Fo,Ml_M0,newman,Mtot_M0\n")
            for r in res:
                f.write("%.8e,%.8e,%.8e,%.10e,%.10e,%.10e\n" % r)


if __name__ == "__main__":
    main()
