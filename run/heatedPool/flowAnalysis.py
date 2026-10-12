#!/usr/bin/env python3
"""
flowAnalysis.py - natural-circulation characterisation of a heated-pool case (latest time).

    python3 flowAnalysis.py <case>          (OpenFOAM environment sourced; parallel cases are reconstructed)

Samples with OpenFOAM (system/flowSampling): a meridional cut through the axis (y = 0, both halves:
theta = 0 at x > 0 and theta = 180 deg at x < 0) and horizontal cuts at z/H = 0.25, 0.5, 0.75. Then writes to
<case>/analysis/:
    meridional.png   T with circulation streamlines and speed with flow-direction arrows (r-z plane)
    radial.png       U_z and T across the liquid at several heights
    axial.png        T along z at five radii (+ thermocouples) and strongest up- and down-flow
    sections.png     T and U_z on horizontal cuts (wedge: radial profiles revolved through 360 deg)
    flowSummary.txt  velocities, circulating mass flow, Re, Pe, Ra, boundary layers, swirl and asymmetry
Velocities in cylindrical components: U_r (outward), U_theta (counter-clockwise), U_z (up). The stream
function assumes axisymmetry: psi(r, z) = rho * integral 2 pi r U_z dr  [kg/s] (mass flow through the
annulus/disc from the inner boundary to r).
"""
import glob
import math
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

FRACS = (0.25, 0.5, 0.75)                       # horizontal cuts, z/H
PROFILE_FRACS = (0.1, 0.25, 0.5, 0.75, 0.9)     # radial profiles, z/H

# light-mode palette (dataviz reference instance)
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
BLUE_RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
WARM_RAMP = ["#fde3d6", "#f9c0a5", "#f39871", "#eb6834", "#c9501f", "#9c3c14", "#6e2a0c"]
DIVERGING = ["#0d366b", "#256abf", "#6da7ec", "#cde2fb", "#f0efec", "#f8d2d1", "#ee8e8d", "#e34948", "#a32726"]
CMAP = "turbo"          # rainbow colormap for all contours (user preference); e.g. "jet" for the classic one


# ----------------------------------------------------------------------------- case data
def load_settings(case):
    import importlib.util
    spec = importlib.util.spec_from_file_location("caseSettings", os.path.join(case, "caseSettings.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def geometry(cfg):
    if cfg.GEOMETRY == "A":
        return dict(ri=cfg.A_INNER_DIAMETER / 2, ro=cfg.A_OUTER_DIAMETER / 2, H=cfg.A_HEIGHT)
    return dict(ri=0.0, ro=cfg.B_DIAMETER / 2, H=cfg.B_HEIGHT)


def latest_time(case):
    times = []
    for d in os.listdir(case):
        try:
            if float(d) > 0 and os.path.isfile(os.path.join(case, d, "T")):
                times.append(d)
        except ValueError:
            pass
    return max(times, key=float) if times else None


def run(cmd, case, log):
    with open(os.path.join(case, log), "w") as fh:
        r = subprocess.run(cmd, cwd=case, stdout=fh, stderr=subprocess.STDOUT)
    if r.returncode:
        raise SystemExit(f"{' '.join(cmd)} failed: see {log}")


def sample(case, g):
    """Reconstruct the latest time if needed and run the cutting-plane sampling."""
    if glob.glob(os.path.join(case, "processor0")):
        run(["reconstructPar", "-latestTime"], case, "log.reconstructPar.flow")
    t = latest_time(case)
    if t is None:
        raise SystemExit("no solution time with T in the case")
    planes = "\n".join(
        f"            z{int(100 * f)} {{ type cuttingPlane; point (0 0 {f * g['H']:.8g}); normal (0 0 1); "
        f"interpolate false; }}" for f in FRACS)
    with open(os.path.join(case, "system", "flowSampling"), "w") as fh:
        fh.write(f"""FoamFile {{ version 2.0; format ascii; class dictionary; object flowSampling; }}
functions
{{
    flowSurfaces
    {{
        type                surfaces;
        libs                (sampling);
        surfaceFormat       raw;
        fields              (T U);
        interpolationScheme cell;
        surfaces
        {{
            meridional {{ type cuttingPlane; point (0 0 0); normal (0 1 0); interpolate false; }}
{planes}
        }}
    }}
}}
""")
    run(["postProcess", "-time", t, "-dict", "system/flowSampling", "-fields", "(T U)"], case, "log.flowSampling")
    return t


def read_raw(case, t, surface, field):
    """Face-centred raw surface data: x y z value(s)."""
    pats = [os.path.join(case, "postProcessing", "flowSurfaces", t, f"{field}_{surface}.raw"),
            os.path.join(case, "postProcessing", "flowSurfaces", t, surface, f"{field}.raw"),
            os.path.join(case, "postProcessing", "flowSurfaces", t, f"{surface}_{field}.raw")]
    for p in pats:
        if os.path.isfile(p):
            return np.loadtxt(p, comments="#")
    found = glob.glob(os.path.join(case, "postProcessing", "flowSurfaces", t, "**", "*"), recursive=True)
    raise SystemExit(f"no raw file for {field} on {surface}; found: {found[:8]}")


# ----------------------------------------------------------------------------- fields on a grid
def half_plane(case, t, side):
    """Meridional half (side = +1: theta = 0, x > 0; -1: theta = 180 deg) as (r, z, T, Ur, Ut, Uz)."""
    T = read_raw(case, t, "meridional", "T")
    U = read_raw(case, t, "meridional", "U")
    m = side * T[:, 0] > 0
    r, z = np.abs(T[m, 0]), T[m, 2]
    s = np.sign(U[m, 0])
    return dict(r=r, z=z, T=T[m, 3], Ur=U[m, 3] * s, Ut=U[m, 4] * s, Uz=U[m, 5])


def to_grid(h, g, nr=60, nz=240):
    """Interpolate scattered face data onto a regular (r, z) grid (linear, on a Delaunay triangulation)."""
    import matplotlib.tri as mtri
    tri = mtri.Triangulation(h["r"], h["z"])
    rr = np.linspace(g["ri"], g["ro"], nr)
    zz = np.linspace(0.0, g["H"], nz)
    R, Z = np.meshgrid(rr, zz)
    out = dict(r=rr, z=zz, R=R, Z=Z)
    for k in ("T", "Ur", "Ut", "Uz"):
        f = mtri.LinearTriInterpolator(tri, h[k])(R, Z)
        a = np.ma.filled(f, np.nan)
        for j in range(a.shape[1]):          # fill the half-cell strips at the walls from the nearest value
            col = a[:, j]
            ok = ~np.isnan(col)
            if ok.any() and (~ok).any():
                a[:, j] = np.interp(zz, zz[ok], col[ok])
        for i in range(a.shape[0]):
            row = a[i, :]
            ok = ~np.isnan(row)
            if ok.any() and (~ok).any():
                a[i, :] = np.interp(rr, rr[ok], row[ok])
        out[k] = a
    return out


def stream_function(gr, rho):
    """psi(r, z) = rho * int_{r_in}^{r} 2 pi r' U_z dr'  [kg/s], axisymmetric."""
    r = gr["r"]
    f = 2 * np.pi * r[None, :] * gr["Uz"]
    psi = np.zeros_like(f)
    psi[:, 1:] = np.cumsum(0.5 * (f[:, 1:] + f[:, :-1]) * np.diff(r)[None, :], axis=1)
    return rho * psi


# ----------------------------------------------------------------------------- plots
def style(ax):
    ax.set_facecolor(SURF)
    ax.grid(True, color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def cmap(colors, name):
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list(name, colors)


def plot_meridional(gr, psi, g, cfg, out):
    """Two panels: temperature, and speed with streamlines carrying flow-direction arrows."""
    import matplotlib.pyplot as plt
    rmm, zmm = 1e3 * gr["R"], 1e3 * gr["Z"]
    TC = gr["T"] - 273.15
    speed = 1e3 * np.sqrt(gr["Ur"] ** 2 + gr["Uz"] ** 2 + gr["Ut"] ** 2)
    fig, axs = plt.subplots(1, 2, figsize=(9, 8), sharey=True, facecolor=SURF)
    cs = axs[0].contourf(rmm, zmm, TC, levels=32, cmap=CMAP)
    fig.colorbar(cs, ax=axs[0], orientation="horizontal", pad=0.07, fraction=0.04).ax.tick_params(labelsize=7)
    axs[0].set_title("Temperature [°C]", fontsize=9, loc="left")
    cs = axs[1].contourf(rmm, zmm, speed, levels=32, cmap=CMAP)
    fig.colorbar(cs, ax=axs[1], orientation="horizontal", pad=0.07, fraction=0.04).ax.tick_params(labelsize=7)
    axs[1].set_title("Speed [mm/s] + streamlines (arrows: flow direction)", fontsize=9, loc="left")
    lev = np.linspace(np.nanmin(psi), np.nanmax(psi), 12)[1:-1]
    lines = axs[1].contour(rmm, zmm, psi, levels=lev, colors="white", linewidths=0.8)
    # arrows ON the streamlines: tangent of the drawn line in display space (the radial axis is stretched),
    # oriented with the local velocity, a few per line at equal display arc length
    for ax in axs:
        ax.set_xlabel("r [mm]")
        style(ax)
    axs[0].set_ylabel("z [mm]")
    fig.suptitle(f"{cfg.GEOMETRY}: meridional plane (radial axis stretched); streamlines = equal mass-flow steps",
                 x=0.01, ha="left", fontsize=10, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.canvas.draw()
    to_disp = axs[1].transData.transform
    rr, zz = 1e3 * gr["r"], 1e3 * gr["z"]
    for segs in lines.allsegs:
        for seg in segs:
            if len(seg) < 8:
                continue
            d = to_disp(seg)
            arc = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(d, axis=0).T))])
            if arc[-1] < 60:                                  # px: too short to carry an arrow
                continue
            n = max(2, int(arc[-1] // 140))
            for sk in (np.arange(n) + 0.5) * arc[-1] / n:
                i = min(int(np.searchsorted(arc, sk)), len(seg) - 2)
                tx, ty = d[i + 1] - d[i]
                if tx == 0 and ty == 0:
                    continue
                ir = int(np.clip(np.searchsorted(rr, seg[i, 0]), 0, len(rr) - 1))
                iz = int(np.clip(np.searchsorted(zz, seg[i, 1]), 0, len(zz) - 1))
                if tx * gr["Ur"][iz, ir] * 1e3 + ty * gr["Uz"][iz, ir] * 1e3 < 0:
                    tx, ty = -tx, -ty
                xy = seg[i]
                tip = axs[1].transData.inverted().transform(d[i] + 9.0 * np.array([tx, ty]) / np.hypot(tx, ty))
                axs[1].annotate("", xy=tip, xytext=xy, arrowprops=dict(arrowstyle="-|>", color="white", lw=0.8,
                                                                         mutation_scale=11, shrinkA=0, shrinkB=0))
    fig.savefig(out, dpi=150, facecolor=SURF)
    plt.close(fig)


def plot_radial(g0, g, out):
    import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4.6), facecolor=SURF)
    for f, c in zip(PROFILE_FRACS, SERIES):
        j = int(round(f * (len(g0["z"]) - 1)))
        a1.plot(1e3 * g0["r"], 1e3 * g0["Uz"][j], color=c, lw=2, label=f"z/H = {f:g}")
        a2.plot(1e3 * g0["r"], g0["T"][j] - 273.15, color=c, lw=2, label=f"z/H = {f:g}")
    a1.axhline(0, color=INK2, lw=0.8)
    a1.set_xlabel("r [mm]")
    a1.set_ylabel("U_z [mm/s]  (up > 0)")
    a1.set_title("Axial velocity across the liquid", loc="left", fontsize=10)
    a2.set_xlabel("r [mm]")
    a2.set_ylabel("T [°C]")
    a2.set_title("Temperature across the liquid", loc="left", fontsize=10)
    a2.legend(frameon=False, fontsize=8)
    for a in (a1, a2):
        style(a)
    fig.tight_layout()
    fig.savefig(out, dpi=150, facecolor=SURF)
    plt.close(fig)


def plot_axial(gr, g, cfg, tc, out):
    import matplotlib.pyplot as plt
    z = 1e3 * gr["z"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 5.4), sharey=True, facecolor=SURF)
    inner = "inner wall" if cfg.GEOMETRY == "A" else "axis"
    for f, c in zip((0.0, 0.25, 0.5, 0.75, 1.0), SERIES):
        j = int(round(f * (len(gr["r"]) - 1)))
        lbl = inner if f == 0 else ("outer wall" if f == 1 else f"r = {1e3 * gr['r'][j]:.2f} mm")
        a1.plot(gr["T"][:, j] - 273.15, z, color=c, lw=2, label=lbl)
    if tc is not None:
        a1.scatter(np.array(tc[1]) - 273.15, 1e3 * np.array(tc[0]), s=40, color=INK, zorder=5,
                   edgecolor=SURF, linewidth=1.5, label="axis thermocouples")
    a1.set_xlabel("T [°C]")
    a1.set_ylabel("z [mm]")
    a1.set_title("Temperature along z at five radii", loc="left", fontsize=10)
    a1.legend(frameon=False, fontsize=8)
    a2.plot(1e3 * np.nanmax(gr["Uz"], axis=1), z, color=SERIES[1], lw=2, label="strongest up-flow")
    a2.plot(1e3 * np.nanmin(gr["Uz"], axis=1), z, color=SERIES[0], lw=2, label="strongest down-flow")
    a2.axvline(0, color=INK2, lw=0.8)
    a2.set_xlabel("U_z [mm/s]")
    a2.set_title("Up- and down-flow along z", loc="left", fontsize=10)
    a2.legend(frameon=False, fontsize=8)
    for a in (a1, a2):
        style(a)
    fig.tight_layout()
    fig.savefig(out, dpi=150, facecolor=SURF)
    plt.close(fig)


def sections_figure(panels, title, out):
    """2 x len(FRACS) section contours, one colour range and one colourbar per row (T, U_z), so the
    heights compare at a glance. panels[row][k] = (X, Y, data, tri or None)."""
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    fig, axs = plt.subplots(2, len(FRACS), figsize=(4.0 * len(FRACS) + 1.2, 8), facecolor=SURF, layout="constrained")
    names = ("T [°C]", "U_z [mm/s]")
    for row in range(2):
        lo = min(float(np.nanmin(pn[2])) for pn in panels[row])
        hi = max(float(np.nanmax(pn[2])) for pn in panels[row])
        pad = 1e-6 * (hi - lo)                  # the extreme value itself must fall inside the top level
        levels = np.linspace(lo - pad, hi + pad, 33)
        for k, f in enumerate(FRACS):
            X, Y, data, tri = panels[row][k]
            ax = axs[row, k]
            if tri is not None:
                cs = ax.tricontourf(tri, data, levels=levels, cmap=CMAP, extend="neither")
            else:
                cs = ax.contourf(X, Y, data, levels=levels, cmap=CMAP)
            ax.set_title(f"{names[row]}, z/H = {f:g}", loc="left", fontsize=9)
            ax.set_aspect("equal")
            ax.set_xlabel("x [mm]")
            ax.set_ylabel("y [mm]")
            style(ax)
        cb = fig.colorbar(cs, ax=list(axs[row, :]), fraction=0.03, pad=0.02,
                          ticks=MaxNLocator(nbins=8).tick_values(lo, hi))
        cb.set_label(f"{names[row]} (same range for all heights)", fontsize=8)
        cb.ax.tick_params(labelsize=7)
    if title:
        fig.suptitle(title, x=0.01, ha="left", fontsize=10, color=INK)
    fig.savefig(out, dpi=140, facecolor=SURF)
    plt.close(fig)


def plot_sections_revolved(gr, g, cfg, out):
    """Wedge (axisymmetric): the radial profiles of T and U_z at each cut height revolved through 360 deg,
    i.e. the wedge repeated all around, in the same layout as the 3D sections."""
    th = np.linspace(0.0, 2 * np.pi, 181)
    R, TH = np.meshgrid(gr["r"], th)
    X, Y = 1e3 * R * np.cos(TH), 1e3 * R * np.sin(TH)
    panels = [[], []]
    for f in FRACS:
        j = int(round(f * (len(gr["z"]) - 1)))
        panels[0].append((X, Y, np.tile(gr["T"][j] - 273.15, (len(th), 1)), None))
        panels[1].append((X, Y, np.tile(1e3 * gr["Uz"][j], (len(th), 1)), None))
    sections_figure(panels, f"{cfg.GEOMETRY} wedge ({cfg.WEDGE_ANGLE:g}° of 360°): radial profiles revolved "
                            "around the axis", out)


def plot_sections(case, t, g, cfg, out):
    """Horizontal cuts: T and U_z contours; the variation around theta goes to the summary only."""
    import matplotlib.tri as mtri
    panels = [[], []]
    stats = []
    for f in FRACS:
        name = f"z{int(100 * f)}"
        T = read_raw(case, t, name, "T")
        U = read_raw(case, t, name, "U")
        x, y = T[:, 0], T[:, 1]
        r = np.hypot(x, y)
        Ur = (U[:, 3] * x + U[:, 4] * y) / np.maximum(r, 1e-12)
        Ut = (-U[:, 3] * y + U[:, 4] * x) / np.maximum(r, 1e-12)
        Uz = U[:, 5]
        tri = mtri.Triangulation(1e3 * x, 1e3 * y)
        if g["ri"] > 0:
            xc = 1e-3 * tri.x[tri.triangles].mean(axis=1)
            yc = 1e-3 * tri.y[tri.triangles].mean(axis=1)
            tri.set_mask(np.hypot(xc, yc) < g["ri"])
        panels[0].append((None, None, T[:, 3] - 273.15, tri))
        panels[1].append((None, None, 1e3 * Uz, tri))
        rm = 0.5 * (g["ri"] + g["ro"])
        thc = np.linspace(-np.pi, np.pi, 361)
        xc, yc = 1e3 * rm * np.cos(thc), 1e3 * rm * np.sin(thc)
        Tc = np.ma.filled(mtri.LinearTriInterpolator(tri, T[:, 3])(xc, yc), np.nan)
        Uc = np.ma.filled(mtri.LinearTriInterpolator(tri, Uz)(xc, yc), np.nan)
        speed = np.sqrt(Ur ** 2 + Ut ** 2 + Uz ** 2)
        stats.append(dict(f=f, T_spread=float(np.nanmax(Tc) - np.nanmin(Tc)),
                          Uz_spread=float(np.nanmax(Uc) - np.nanmin(Uc)),
                          Ut_max=float(np.max(np.abs(Ut))), U_max=float(np.max(speed)),
                          Ur_max=float(np.max(np.abs(Ur)))))
    sections_figure(panels, None, out)
    return stats


# ----------------------------------------------------------------------------- summary
def summary(cfg, g, g0, g180, psi, stats, t, case):
    rho, mu, k, cp, beta = cfg.RHO, None, cfg.K, cfg.CP, cfg.BETA
    vis = cfg.VISCOSITY
    mu = vis["mu0"] * math.exp(vis["Ea"] / (8.314462618 * cfg.T_REF)) if vis["model"] == "arrhenius" else cfg.MU
    nu, alpha = mu / rho, k / (rho * cp)
    L = g["ro"] - g["ri"] if cfg.GEOMETRY == "A" else 2 * g["ro"]
    speed = np.sqrt(g0["Ur"] ** 2 + g0["Ut"] ** 2 + g0["Uz"] ** 2)
    Umax = float(np.nanmax(speed))
    up, down = float(np.nanmax(g0["Uz"])), float(np.nanmin(g0["Uz"]))
    circ = float(np.nanmax(np.nanmax(psi, axis=1) - np.nanmin(psi, axis=1)) / 2)
    j = len(g0["z"]) // 2
    prof = g0["Uz"][j]
    # down-flow along the cooled outer wall: distance from the wall to its strongest point, and the width
    # of the down-flow layer (wall to where U_z changes sign)
    i_dn = int(np.nanargmin(prof))
    d_vel = g["ro"] - g0["r"][i_dn]
    up_side = np.where(prof[:i_dn] > 0)[0]
    d_layer = g["ro"] - (g0["r"][up_side[-1]] if up_side.size else g["ri"])
    dT_rad = float(np.nanmax(g0["T"][j]) - np.nanmin(g0["T"][j]))
    dT_ax = float(np.nanmax(g0["T"]) - np.nanmin(g0["T"]))
    Ra = 9.81 * beta * dT_ax * g["H"] ** 3 / (nu * alpha)
    asym = None
    if g180 is not None:
        inner = (slice(3, -3), slice(3, -3))   # interior: wall strips are interpolation-filled
        # percentiles: single grid points in the steep wall layers carry interpolation noise
        asym = (float(np.nanpercentile(np.abs(g0["T"][inner] - g180["T"][inner]), 90)),
                float(np.nanpercentile(np.abs(g0["Uz"][inner] - g180["Uz"][inner]), 90)))
    lines = [
        f"Natural circulation, case {os.path.basename(os.path.abspath(case))}, time {t}",
        f"  liquid: nu = {nu:.3g} m2/s, alpha = {alpha:.3g} m2/s, Pr = {nu / alpha:.2f}; length L = {1e3 * L:.2f} mm "
        f"({'gap' if cfg.GEOMETRY == 'A' else 'diameter'}), H = {1e3 * g['H']:.1f} mm",
        f"  max speed            : {1e3 * Umax:.6g} mm/s",
        f"  strongest up / down  : {1e3 * up:.6g} / {1e3 * down:.6g} mm/s",
        f"  circulating mass flow: {1e6 * circ:.6g} mg/s  (each way; axisymmetric stream function)",
        f"  turnover time H/U_max: {g['H'] / max(Umax, 1e-12):.3g} s",
        f"  Re = U_max L / nu    : {Umax * L / nu:.3g}",
        f"  Pe = U_max L / alpha : {Umax * L / alpha:.3g}  (> 1: heat carried by the flow, < 1: conduction)",
        f"  Ra_H (dT over height): {Ra:.3g}  (dT = {dT_ax:.1f} K, H)",
        f"  dT across the liquid at mid-height: {dT_rad:.4f} K",
        f"  down-flow at the outer wall (mid-height): strongest {1e3 * d_vel:.2f} mm from the wall, "
        f"layer {1e3 * d_layer:.2f} mm wide (of {1e3 * (g['ro'] - g['ri']):.2f} mm)",
    ]
    if asym:
        lines.append(f"  theta = 0 vs 180 deg (meridional halves, interior, 90th percentile): |dT| = {asym[0]:.2f} K, "
                     f"|dU_z| = {1e3 * asym[1]:.3g} mm/s")
    else:
        lines.append("  axisymmetric wedge: no theta information")
    for s in stats:
        lines.append(f"  z/H = {s['f']:.2f}: around theta at mid-radius dT = {s['T_spread']:.2f} K, "
                     f"dU_z = {1e3 * s['Uz_spread']:.3g} mm/s; max |U_theta| = {1e3 * s['Ut_max']:.3g} mm/s "
                     f"({100 * s['Ut_max'] / max(s['U_max'], 1e-12):.1f} % of max speed), max |U_r| = {1e3 * s['Ur_max']:.3g} mm/s")
    return lines


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    case = os.path.abspath(sys.argv[1])
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": GRID, "axes.labelcolor": INK2,
                         "xtick.color": INK2, "ytick.color": INK2, "text.color": INK})
    cfg = load_settings(case)
    g = geometry(cfg)
    t = sample(case, g)
    out = os.path.join(case, "analysis")
    os.makedirs(out, exist_ok=True)
    wedge = getattr(cfg, "MESH_TYPE", "3D") == "wedge"
    g0 = to_grid(half_plane(case, t, +1), g)
    g180 = None if wedge else to_grid(half_plane(case, t, -1), g)
    psi = stream_function(g0, cfg.RHO)
    tc = None
    data = (getattr(cfg, "THERMOCOUPLES", {}) or {}).get(cfg.GEOMETRY) or {}
    if "axis" in data:
        tc = (data["axis"]["z"], data["axis"]["T"])
    plot_meridional(g0, psi, g, cfg, os.path.join(out, "meridional.png"))
    plot_radial(g0, g, os.path.join(out, "radial.png"))
    plot_axial(g0, g, cfg, tc, os.path.join(out, "axial.png"))
    figs = "meridional.png, radial.png, axial.png, sections.png"
    stats = []
    if wedge:                         # axisymmetric: revolve the radial profiles (the wedge repeated all around)
        plot_sections_revolved(g0, g, cfg, os.path.join(out, "sections.png"))
    else:
        stats = plot_sections(case, t, g, cfg, os.path.join(out, "sections.png"))
    lines = summary(cfg, g, g0, g180, psi, stats, t, case)
    with open(os.path.join(out, "flowSummary.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nwrote {out}/{figs}, flowSummary.txt")


if __name__ == "__main__":
    main()
