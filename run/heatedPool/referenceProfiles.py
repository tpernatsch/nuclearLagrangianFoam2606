#!/usr/bin/env python3
"""
referenceProfiles.py - shapes and heat contributions from a reference model of a similar capsule
(SALIENT-03 "Salt surface temp and flux.xlsx": salt-side wall temperature and heat flux on the side
wall and the bottom, steady state).

    python3 referenceProfiles.py [xlsx]     tables on screen, analysis/reference_profiles.png,
                                            analysis/reference_shapes.csv (z/H, theta, q/q_mean)

Columns: x = axial (our z), y = radial. The side wall is listed at y = 0.038 but the bottom ends at
y = 0.0038 and its edge temperature equals the side's bottom temperature, so the radius is 3.8 mm
(units slip in the side's y). Only shapes and contributions are used, not absolute values.
"""
import csv
import math
import os
import re
import sys
import xml.etree.ElementTree as ET
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
XLSX = os.path.expanduser("/mnt/c/Users/tomma/Downloads/PhD/PeriodAbroad/SALIENT03/Salt surface temp and flux.xlsx")
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
R_REF = 3.8e-3                       # reference capsule radius [m] (see above)

# our A axis thermocouples, for a shape comparison (z [m], T [K])
TC_A = ([0.0, 37.5e-3, 75.0e-3], [893.15, 963.15, 983.15])


def read_sheet(path):
    """Cells of the first sheet as {row: {col: value}} (standard library only)."""
    with zipfile.ZipFile(path) as z:
        ss = []
        if "xl/sharedStrings.xml" in z.namelist():
            for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", NS):
                ss.append("".join(si.itertext()))
        root = ET.fromstring(z.read("xl/worksheets/sheet1.xml"))
    rows = {}
    for c in root.iter("{%s}c" % NS["m"]):
        v = c.find("m:v", NS)
        if v is None:
            continue
        col, row = re.match(r"([A-Z]+)(\d+)", c.get("r")).groups()
        rows.setdefault(int(row), {})[col] = ss[int(v.text)] if c.get("t") == "s" else v.text
    return rows


def series(rows, cx, cy, cv):
    out = []
    for r in sorted(rows):
        d = rows[r]
        if r > 1 and cx in d and cv in d:
            out.append((float(d[cx]), float(d[cy]), float(d[cv])))
    return out


def trap(xs, fs):
    return sum(0.5 * (fs[i] + fs[i + 1]) * (xs[i + 1] - xs[i]) for i in range(len(xs) - 1))


def interp(xs, ys, x):
    if x <= xs[0]:
        return ys[0]
    for i in range(1, len(xs)):
        if x <= xs[i]:
            return ys[i - 1] + (ys[i] - ys[i - 1]) * (x - xs[i - 1]) / (xs[i] - xs[i - 1])
    return ys[-1]


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else XLSX
    rows = read_sheet(path)
    ts, tb = series(rows, "A", "B", "C"), series(rows, "E", "F", "G")
    qs, qb = series(rows, "I", "J", "K"), series(rows, "M", "N", "O")

    zt, Tt = [p[0] for p in ts], [p[2] for p in ts]
    zq, qq = [p[0] for p in qs], [p[2] for p in qs]
    rb, Tb = [p[1] for p in tb], [p[2] for p in tb]
    rq, qbv = [p[1] for p in qb], [p[2] for p in qb]
    H = max(zt[-1], zq[-1])
    Rb = rb[-1]

    Qs = 2 * math.pi * R_REF * trap(zq, qq)
    Qb = 2 * math.pi * trap(rq, [r * q for r, q in zip(rq, qbv)])
    V = math.pi * R_REF ** 2 * H
    print(f"Reference capsule: R = {1e3 * R_REF:.1f} mm, H = {1e3 * H:.1f} mm (bottom data to r = {1e3 * Rb:.2f} mm)")
    print(f"  side   : {Qs:7.1f} W  ({100 * Qs / (Qs + Qb):.1f} %), mean {Qs / (2 * math.pi * R_REF * H) / 1e3:.0f} kW/m2")
    print(f"  bottom : {Qb:7.1f} W  ({100 * Qb / (Qs + Qb):.1f} %), mean {Qb / (math.pi * Rb ** 2) / 1e3:.0f} kW/m2")
    print(f"  total  : {Qs + Qb:7.1f} W  -> {(Qs + Qb) / V / 1e6:.0f} W/cm3 if it is all salt generation (top adiabatic)")
    print(f"  wall T : side {min(Tt):.1f} .. {max(Tt):.1f} K (max at z/H = {zt[Tt.index(max(Tt))] / H:.2f}), "
          f"bottom {min(Tb):.1f} .. {max(Tb):.1f} K")

    # cumulative share of the side heat from the bottom up
    cum = [0.0]
    for i in range(1, len(zq)):
        cum.append(cum[-1] + 0.5 * (qq[i] + qq[i - 1]) * (zq[i] - zq[i - 1]))
    print("\n  side heat released below z/H:  " + ",  ".join(
        f"{f:.2f}: {100 * interp(zq, cum, f * H) / cum[-1]:.0f} %" for f in (0.25, 0.5, 0.75, 0.9, 0.95)))
    print(f"  top 10 % of the height releases {100 * (1 - interp(zq, cum, 0.9 * H) / cum[-1]):.0f} % of the side heat "
          f"(flux at the top edge {qq[-1] / 1e3:.0f} kW/m2)")

    # normalised shapes on a common grid
    Tmin, Tmax = min(Tt), max(Tt)
    qmean = cum[-1] / (zq[-1] - zq[0])
    grid = [i / 100 for i in range(101)]
    shape = [(f, (interp(zt, Tt, f * H) - Tmin) / (Tmax - Tmin), interp(zq, qq, f * H) / qmean) for f in grid]
    print("\n   z/H   theta=(T-Tmin)/(Tmax-Tmin)   q/q_mean")
    for f, th, qn in shape[::10]:
        print(f"  {f:5.2f}   {th:8.3f}                     {qn:6.2f}")
    # local q vs T: is the side loss a single Robin law?
    print("\n  local slope dq/dT along the side [kW/m2/K] (constant would mean one Robin law q = H (T - T_ref)):")
    for a, b in ((0.1, 0.3), (0.3, 0.5), (0.5, 0.7), (0.7, 0.85)):
        dq = interp(zq, qq, b * H) - interp(zq, qq, a * H)
        dT = interp(zt, Tt, b * H) - interp(zt, Tt, a * H)
        print(f"    z/H {a:.2f}-{b:.2f}: {dq / dT / 1e3:6.2f}")

    os.makedirs(os.path.join(HERE, "analysis"), exist_ok=True)
    with open(os.path.join(HERE, "analysis", "reference_shapes.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["z_over_H", "theta_side", "q_side_over_mean"])
        for f, th, qn in shape:
            w.writerow([f"{f:.2f}", f"{th:.5f}", f"{qn:.5f}"])
    plot(H, zt, Tt, zq, qq, qmean, rb, Tb, rq, qbv, cum)


def plot(H, zt, Tt, zq, qq, qmean, rb, Tb, rq, qbv, cum):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    blue, orange, aqua = "#2a78d6", "#eb6834", "#1baf7a"     # categorical slots 1-3 (light)
    ink, ink2, grid, surf = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": grid, "axes.labelcolor": ink2,
                         "xtick.color": ink2, "ytick.color": ink2, "text.color": ink})
    fig, axs = plt.subplots(1, 3, figsize=(11, 5), sharey=True, facecolor=surf)
    Tmin, Tmax = min(Tt), max(Tt)

    ax = axs[0]
    ax.plot([(t - Tmin) / (Tmax - Tmin) for t in Tt], [z / H for z in zt], color=blue, lw=2,
            label="reference, side wall (salt side)")
    tz, tT = TC_A
    lo, hi = min(tT), max(tT)
    ax.scatter([(t - lo) / (hi - lo) for t in tT], [z / tz[-1] for z in tz], s=46, color=ink, zorder=5,
               edgecolor=surf, linewidth=2, label="our A axis thermocouples")
    ax.set_xlabel("θ = (T − T_min)/(T_max − T_min)")
    ax.set_ylabel("z / H")
    ax.set_title("Wall temperature shape", loc="left", fontsize=10)
    ax.legend(frameon=False, loc="lower right", fontsize=8)

    ax = axs[1]
    ax.plot([q / qmean for q in qq], [z / H for z in zq], color=orange, lw=2)
    ax.axvline(1.0, color=ink2, lw=0.8, ls=(0, (3, 3)))
    ax.annotate("mean", (1.0, 0.03), xytext=(4, 0), textcoords="offset points", color=ink2, fontsize=8)
    ax.set_xlabel("q / q_mean  (side wall)")
    ax.set_title("Side heat flux shape", loc="left", fontsize=10)
    ax.set_xlim(0, 4)

    ax = axs[2]
    ax.plot([c / cum[-1] for c in cum], [z / H for z in zq], color=aqua, lw=2)
    ax.set_xlabel("share of the side heat released below z")
    ax.set_title("Cumulative side heat", loc="left", fontsize=10)
    ax.set_xlim(0, 1)
    for a in axs:
        a.set_facecolor(surf)
        a.grid(True, color=grid, lw=0.6)
        a.set_axisbelow(True)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    axs[0].set_ylim(0, 1.0)
    fig.suptitle("Reference capsule (SALIENT-03 model, R 3.8 mm, H 36 mm): shapes for our boundary conditions",
                 x=0.01, ha="left", fontsize=10, color=ink)
    fig.tight_layout()
    out = os.path.join(HERE, "analysis", "reference_profiles.png")
    fig.savefig(out, dpi=150, facecolor=surf)
    print(f"\nwrote {out} and analysis/reference_shapes.csv")


if __name__ == "__main__":
    main()
