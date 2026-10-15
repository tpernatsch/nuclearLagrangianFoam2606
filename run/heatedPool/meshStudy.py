#!/usr/bin/env python3
"""
meshStudy.py - mesh-sensitivity study on axisymmetric wedges (the steady 3D flows are axisymmetric).

    python3 meshStudy.py run        generate + run + post-process all levels (OpenFOAM sourced)
    python3 meshStudy.py analyse    tables, Grid Convergence Index and plots from finished runs

Cases go to cases/meshStudy/<geometry>_<level>_wedge. Every part of the mesh is refined together: the target cell
sizes of the level (setupCase.MESH_LEVELS), the first boundary-layer cell (30 um at medium) and the layer growth
ratio (1.15 at medium, 1.15^(h/h_medium)), all scaled with the level (30 um at
medium). Quantities of interest come from postProcessing/heatBalance_*.dat and analysis/flowSummary.txt.
GCI after Celik et al. (2008), J. Fluids Eng. 130, 078001, on the three finest levels; the representative
cell size of the 2D (r-z) mesh is h = sqrt(A_rz / N_cells).
Output: cases/meshStudy/meshStudy.txt, meshStudy.csv, meshStudy.png
"""
import csv
import math
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import setupCase as sc  # noqa: E402

LEVELS = ["coarse", "medium", "fine", "xfine", "xxfine"]
GEOS = ["A", "B"]
BASE = os.path.join(HERE, "cases", "meshStudy")
BL_MEDIUM = 30e-6
MAX_PARALLEL = 10

# quantity: (label, unit, where it comes from)
QOIS = {
    "Tavg":    ("mean salt T", "K"),
    "Tmax":    ("max salt T", "K"),
    "Umax":    ("max speed", "mm/s"),
    "Uup":     ("strongest up-flow", "mm/s"),
    "Udown":   ("strongest down-flow", "mm/s"),
    "circ":    ("circulating mass flow", "mg/s"),
    "dTmid":   ("dT across the salt at mid-height", "K"),
    "Qring":   ("ring heat into the salt (A)", "W"),
    "Qouter":  ("heat through the outer wall", "W"),
}


def case_dir(geo, lvl):
    return os.path.join(BASE, f"{geo}_{lvl}_wedge")      # setupCase appends _wedge


def config_text(geo, lvl):
    h = sc.MESH_LEVELS[lvl]["h"] / sc.MESH_LEVELS["medium"]["h"]
    return (f'\nGEOMETRY = "{geo}"\nMESH_LEVEL = "{lvl}"\nMESH_TYPE = "wedge"\nBOUNDARY_LAYER = True\n'
            f'BL_FIRST_CELL = {BL_MEDIUM * h:.6g}\nBL_GROWTH = {1.15 ** h:.6g}\n'
            f'RUN_MODE = "steady"\nSTEADY_ITERATIONS = 40000\n'
            f'STEADY_WRITE_INTERVAL = 40000\nHEAT_BALANCE_EVERY = 100\n'
            f'OUTPUT_DIR = "{BASE}"\nCASE_NAME = "{geo}_{lvl}"\nN_PROCS = 1\n')


def run():
    os.makedirs(BASE, exist_ok=True)
    defaults = open(os.path.join(HERE, "caseDefaults.py")).read()
    for geo in GEOS:
        for lvl in LEVELS:
            cfgfile = os.path.join(BASE, f"settings_{geo}_{lvl}.py")
            with open(cfgfile, "w") as fh:
                fh.write(defaults + config_text(geo, lvl))
            subprocess.run([sys.executable, os.path.join(HERE, "setupCase.py"), "--from", cfgfile, "--force"],
                           check=True, stdout=subprocess.DEVNULL)
    jobs = [case_dir(g, l) for g in GEOS for l in LEVELS]
    running, t0 = [], time.time()
    while jobs or running:
        while jobs and len(running) < MAX_PARALLEL:
            c = jobs.pop(0)
            running.append((c, subprocess.Popen(["./Allrun"], cwd=c, stdout=open(os.path.join(c, "allrun.log"), "w"),
                                                stderr=subprocess.STDOUT)))
        for item in list(running):
            if item[1].poll() is not None:
                running.remove(item)
                print(f"  done {os.path.basename(item[0])} ({time.time() - t0:.0f} s)", flush=True)
        time.sleep(2)
    for c in (case_dir(g, l) for g in GEOS for l in LEVELS):
        subprocess.run([sys.executable, os.path.join(HERE, "flowAnalysis.py"), c], stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
    analyse()


# ----------------------------------------------------------------------------- results
def read_case(geo, lvl):
    c = case_dir(geo, lvl)
    out = dict(geo=geo, level=lvl)
    log = os.path.join(c, "log.buoyantBoussinesqSimpleFoam.steady")
    txt = open(log).read() if os.path.isfile(log) else ""
    m = re.search(r"SIMPLE solution converged in (\d+) iterations", txt)
    out["converged"] = bool(m)
    out["iterations"] = int(m.group(1)) if m else len(re.findall(r"^Time = ", txt, re.M))
    cm = re.search(r"cells:\s+(\d+)", open(os.path.join(c, "log.checkMesh")).read())
    out["cells"] = int(cm.group(1))
    hb = [f for f in os.listdir(os.path.join(c, "postProcessing")) if f.startswith("heatBalance_")][0]
    lines = [l.split() for l in open(os.path.join(c, "postProcessing", hb)) if l.strip()]
    head = lines[0][1:]
    last = dict(zip(head, (float(v) for v in lines[-1])))
    out["Tavg"], out["Tmax"] = last["Tavg"], last["Tmax"]
    out["Qring"] = -last["Q_innerWall"] if "Q_innerWall" in last else None
    out["Qouter"] = last["Q_outerWall"]
    fs = open(os.path.join(c, "analysis", "flowSummary.txt")).read()
    num = r"([-+0-9.eE]+)"
    out["Umax"] = float(re.search(r"max speed\s*:\s*" + num, fs).group(1))
    up, down = re.search(r"strongest up / down\s*:\s*" + num + r"\s*/\s*" + num, fs).groups()
    out["Uup"], out["Udown"] = float(up), -float(down)
    out["circ"] = float(re.search(r"circulating mass flow:\s*" + num, fs).group(1))
    out["dTmid"] = float(re.search(r"dT across the liquid at mid-height:\s*" + num, fs).group(1))
    cfg = sc.load_config(os.path.join(c, "caseSettings.py"))
    g = (cfg.A_OUTER_DIAMETER - cfg.A_INNER_DIAMETER) / 2 if geo == "A" else cfg.B_DIAMETER / 2
    H = cfg.A_HEIGHT if geo == "A" else cfg.B_HEIGHT
    out["h"] = math.sqrt(g * H / out["cells"])          # representative cell size of the r-z mesh
    return out


def gci(phi, h):
    """Celik et al. (2008): phi, h ordered fine -> coarse (3 levels). Returns p, phi_ext, e_a, GCI_fine, note."""
    f1, f2, f3 = phi
    r21, r32 = h[1] / h[0], h[2] / h[1]
    e21, e32 = f2 - f1, f3 - f2
    if abs(e21) < 1e-14 * max(abs(f1), 1) or abs(e32) < 1e-14 * max(abs(f1), 1):
        return None, f1, 0.0, 0.0, "converged to round-off"
    s = 1.0 if e32 / e21 > 0 else -1.0
    p = 2.0
    for _ in range(200):
        q = math.log((r21 ** p - s) / (r32 ** p - s)) if (r21 ** p - s) > 0 and (r32 ** p - s) > 0 else 0.0
        p_new = abs(math.log(abs(e32 / e21)) + q) / math.log(r21)
        if abs(p_new - p) < 1e-8:
            p = p_new
            break
        p = p_new
    ext = (r21 ** p * f1 - f2) / (r21 ** p - 1)
    ea = abs((f1 - f2) / f1)
    g_fine = 1.25 * ea / (r21 ** p - 1)
    note = "monotonic" if s > 0 else "oscillatory"
    return p, ext, ea, g_fine, note


def analyse():
    rows = []
    for geo in GEOS:
        for lvl in LEVELS:
            try:
                rows.append(read_case(geo, lvl))
            except Exception as exc:                # noqa: BLE001
                print(f"  {geo} {lvl}: not available ({exc})")
    with open(os.path.join(BASE, "meshStudy.csv"), "w", newline="") as fh:
        keys = ["geo", "level", "cells", "h", "iterations", "converged"] + list(QOIS)
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    lines = ["Mesh-sensitivity study (axisymmetric wedges, steady); GCI on the three finest levels (Celik 2008)", ""]
    for geo in GEOS:
        rs = [r for r in rows if r["geo"] == geo]
        if not rs:
            continue
        lines.append(f"Geometry {geo}")
        hdr = f"  {'level':<8}{'cells':>8}{'h [um]':>8}{'iter':>7}  " + "".join(f"{k:>10}" for k in QOIS)
        lines.append(hdr)
        for r in rs:
            vals = "".join(f"{(r[k] if r[k] is not None else float('nan')):>10.4g}" for k in QOIS)
            lines.append(f"  {r['level']:<8}{r['cells']:>8}{1e6 * r['h']:>8.1f}{r['iterations']:>7}"
                         f"{'' if r['converged'] else '*'}  {vals}")
        fin = sorted(rs, key=lambda r: r["h"])[:3]
        if len(fin) == 3:
            lines.append(f"  GCI ({', '.join(r['level'] for r in fin)}):")
            for k, (lab, unit) in QOIS.items():
                if fin[0][k] is None:
                    continue
                p, ext, ea, g, note = gci([r[k] for r in fin], [r["h"] for r in fin])
                ptxt = f"p = {p:4.2f}" if p is not None else "p =  -  "
                lines.append(f"    {lab:<34} finest {fin[0][k]:10.5g} {unit:<5} {ptxt}  extrapolated {ext:10.5g}"
                             f"  change vs next {100 * ea:6.3f} %  GCI_fine {100 * g:6.3f} %  ({note})")
        lines.append("")
    lines.append("* = stopped at the iteration cap without meeting the residual targets")
    txt = "\n".join(lines)
    with open(os.path.join(BASE, "meshStudy.txt"), "w") as fh:
        fh.write(txt + "\n")
    print(txt)
    plot(rows)


def plot(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ink2, grid, surf = "#52514e", "#e4e3df", "#fcfcfb"
    colors = {"A": "#2a78d6", "B": "#eb6834"}
    keys = [k for k in QOIS]
    n = len(keys)
    cols = 3
    fig, axs = plt.subplots(math.ceil(n / cols), cols, figsize=(12, 3.2 * math.ceil(n / cols)), facecolor=surf)
    for ax, k in zip(axs.flat, keys):
        lab, unit = QOIS[k]
        for geo in GEOS:
            rs = sorted([r for r in rows if r["geo"] == geo and r[k] is not None], key=lambda r: r["h"])
            if not rs:
                continue
            ref = rs[0][k]
            x = [1e6 * r["h"] for r in rs]
            y = [100 * (r[k] - ref) / ref for r in rs]
            ax.plot(x, y, "o-", color=colors[geo], lw=2, ms=6, label=f"{geo} (finest {ref:.4g} {unit})")
        ax.axhline(0, color=ink2, lw=0.7)
        ax.set_title(lab, loc="left", fontsize=9)
        ax.set_xlabel("representative cell size h [um]")
        ax.set_ylabel("difference to finest [%]")
        ax.legend(frameon=False, fontsize=7.5)
        ax.set_facecolor(surf)
        ax.grid(True, color=grid, lw=0.5)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    for ax in list(axs.flat)[n:]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(os.path.join(BASE, "meshStudy.png"), dpi=140, facecolor=surf)
    plt.close(fig)
    print(f"\nwrote {BASE}/meshStudy.txt, meshStudy.csv, meshStudy.png")


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in ("run", "analyse"):
        raise SystemExit(__doc__)
    run() if sys.argv[1] == "run" else analyse()
