#!/usr/bin/env python3
"""
postProcessLDA.py - reduce the run-time outputs of a bubbly-tube case to comparison tables.

    python3 postProcessLDA.py <caseDir> [--tavg T] [--plot]

Writes <caseDir>/postProcessing/summary/:
  planes.csv      per LDA height: area-averaged void (instantaneous mean/std over t >= tAvg and of the
                  mean field), liquid up/down volumetric flow through the plane from the time-averaged
                  velocity, the areas and area-averaged velocities of the up/down regions
                  ([dW] {u}_d analogue), and the same down-flow integrated from the LDA points only
                  (PCHIP along each radial path, [dW] Eq. 6.9)
  lda_points.csv  per LDA point: mean U_z, sigma_z, TI_z = sigma_z/|mean|, mean alpha, samples
  lda_radial.csv  per height and radius: angle-averaged mean U_z, sigma_z and alpha
  monitors.txt    holdup, free-surface level / level swell, gas mass balance
Requires only numpy (matplotlib for --plot).
"""
import glob
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


# ----------------------------------------------------------------------------- readers
def time_dirs(fo_dir):
    out = []
    for d in glob.glob(os.path.join(fo_dir, "*")):
        try:
            out.append((float(os.path.basename(d)), d))
        except ValueError:
            pass
    return [d for _, d in sorted(out)]


def read_dat(fo_dir, fname):
    """Concatenate a function-object .dat file over restart directories -> (t, columns)."""
    rows = {}
    for d in time_dirs(fo_dir):
        p = os.path.join(d, fname)
        if not os.path.exists(p):
            continue
        with open(p) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                v = [float(x) for x in line.split()]
                rows[v[0]] = v[1:]
    if not rows:
        return None, None
    t = np.array(sorted(rows))
    return t, np.array([rows[k] for k in t])


def read_probes(fo_dir, field):
    """Probe file -> (t, data[nt, nprobe, ncomp])."""
    rows = {}
    for d in time_dirs(fo_dir):
        p = os.path.join(d, field)
        if not os.path.exists(p):
            continue
        with open(p) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                t, rest = line.split(None, 1)
                if "(" in rest:
                    vals = [[float(x) for x in g.split()] for g in re.findall(r"\(([^)]*)\)", rest)]
                else:
                    vals = [[float(x)] for x in rest.split()]
                rows[float(t)] = vals
    if not rows:
        return None, None
    t = np.array(sorted(rows))
    return t, np.array([rows[k] for k in t])


def read_config(case):
    import importlib.util
    p = os.path.join(case, "caseSettings.py")
    spec = importlib.util.spec_from_file_location("cfg", p)
    cfg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cfg)
    return cfg


# ----------------------------------------------------------------------------- PCHIP (Fritsch-Carlson)
def pchip(x, y, xi):
    x, y = np.asarray(x, float), np.asarray(y, float)
    h = np.diff(x)
    delta = np.diff(y) / h
    n = len(x)
    m = np.zeros(n)
    if n == 2:
        m[:] = delta[0]
    else:
        for k in range(1, n - 1):
            if delta[k - 1] * delta[k] > 0:
                w1, w2 = 2 * h[k] + h[k - 1], h[k] + 2 * h[k - 1]
                m[k] = (w1 + w2) / (w1 / delta[k - 1] + w2 / delta[k])
        for k, (d0, d1, h0, h1) in ((0, (delta[0], delta[1], h[0], h[1])),
                                     (n - 1, (delta[-1], delta[-2], h[-1], h[-2]))):
            mk = ((2 * h0 + h1) * d0 - h0 * d1) / (h0 + h1)
            if np.sign(mk) != np.sign(d0):
                mk = 0.0
            elif np.sign(d0) != np.sign(d1) and abs(mk) > abs(3 * d0):
                mk = 3 * d0
            m[k] = mk
    xi = np.asarray(xi, float)
    i = np.clip(np.searchsorted(x, xi) - 1, 0, n - 2)
    t = (xi - x[i]) / h[i]
    h00, h10 = 2 * t**3 - 3 * t**2 + 1, t**3 - 2 * t**2 + t
    h01, h11 = -2 * t**3 + 3 * t**2, t**3 - t**2
    return h00 * y[i] + h10 * h[i] * m[i] + h01 * y[i + 1] + h11 * h[i] * m[i + 1]


def lda_flow(radii, angles_deg, W, R):
    """[dW] Eq. 6.9 over a full circle: w(0) from the innermost point, w(R) = 0, PCHIP along each path.
    W[ia, ir] = mean U_z.  Returns (Q_down, Q_up, A_down)."""
    th = np.radians(np.asarray(angles_deg, float))
    order = np.argsort(th)
    th, W = th[order], W[order]
    na = len(th)
    dth = np.empty(na)
    for k in range(na):       # periodic central differences
        dth[k] = 0.5 * ((th[(k + 1) % na] - th[k - 1]) % (2 * np.pi))
    r_fine = np.linspace(0, R, 2001)
    Qd = Qu = Ad = 0.0
    for k in range(na):
        rk = np.concatenate(([0.0], radii, [R]))
        wk = np.concatenate(([W[k, 0]], W[k], [0.0]))
        w = pchip(rk, wk, r_fine)
        Qd += dth[k] * np.trapezoid(np.minimum(w, 0) * r_fine, r_fine)
        Qu += dth[k] * np.trapezoid(np.maximum(w, 0) * r_fine, r_fine)
        Ad += dth[k] * np.trapezoid((w < 0) * r_fine, r_fine)
    return Qd, Qu, Ad


# ----------------------------------------------------------------------------- main
def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if not args:
        raise SystemExit(__doc__)
    case = os.path.abspath(args[0])
    cfg = read_config(case)
    tavg = cfg.AVERAGE_START
    if "--tavg" in sys.argv:
        tavg = float(sys.argv[sys.argv.index("--tavg") + 1])
    pp = os.path.join(case, "postProcessing")
    out = os.path.join(pp, "summary")
    os.makedirs(out, exist_ok=True)
    R = cfg.TUBE_DIAMETER / 2
    A = math.pi * R * R
    Z, AN, RA = cfg.LDA_HEIGHTS, cfg.LDA_ANGLES, cfg.LDA_RADII
    nz, na, nr = len(Z), len(AN), len(RA)

    # ---- probes
    tU, U = read_probes(os.path.join(pp, "ldaProbes"), "U.water")
    ta, al = read_probes(os.path.join(pp, "ldaProbes"), "alpha.air")
    Wmean = np.full((nz, na, nr), np.nan)
    Wstd = np.full((nz, na, nr), np.nan)
    Amean = np.full((nz, na, nr), np.nan)
    nsamp = 0
    if tU is not None and na > 0:
        sel = tU >= tavg
        nsamp = int(sel.sum())
        uz = U[sel, :, 2].reshape(nsamp, nz, na, nr)
        Wmean, Wstd = uz.mean(0), uz.std(0, ddof=1) if nsamp > 1 else np.zeros_like(uz[0])
        if ta is not None:
            Amean = al[ta >= tavg, :, 0].reshape(-1, nz, na, nr).mean(0)
    with open(os.path.join(out, "lda_points.csv"), "w") as f:
        f.write("z,angle_deg,r,x,y,Uz_mean,sigma_z,TI_z,alpha_mean,n_samples\n")
        for iz, z in enumerate(Z):
            for ia, a in enumerate(AN):
                for ir, r in enumerate(RA):
                    m, s = Wmean[iz, ia, ir], Wstd[iz, ia, ir]
                    ti = s / abs(m) if abs(m) > 1e-12 else float("nan")
                    x, y = r * math.cos(math.radians(a)), r * math.sin(math.radians(a))
                    f.write(f"{z},{a},{r},{x:.6g},{y:.6g},{m:.6g},{s:.6g},{ti:.6g},{Amean[iz,ia,ir]:.6g},{nsamp}\n")
    with open(os.path.join(out, "lda_radial.csv"), "w") as f:
        f.write("z,r,Uz_mean,sigma_z,alpha_mean\n")
        for iz, z in enumerate(Z):
            for ir, r in enumerate(RA):
                f.write(f"{z},{r},{np.nanmean(Wmean[iz,:,ir]):.6g},{np.nanmean(Wstd[iz,:,ir]):.6g},"
                        f"{np.nanmean(Amean[iz,:,ir]):.6g}\n")

    # ---- planes
    rows = []
    for iz, z in enumerate(Z, 1):
        t, v = read_dat(os.path.join(pp, f"voidZ{iz}"), "surfaceFieldValue.dat")
        vi_m = vi_s = float("nan")
        if t is not None and (t >= tavg).any():
            vi_m, vi_s = v[t >= tavg, 0].mean(), v[t >= tavg, 0].std()
        t, v = read_dat(os.path.join(pp, f"meanVoidZ{iz}"), "surfaceFieldValue.dat")
        vm = v[-1, 0] if t is not None else float("nan")
        t, v = read_dat(os.path.join(pp, f"flowZ{iz}"), "surfaceFieldValue.dat")
        if t is not None:
            Qd, Qu, Ad = v[-1, 0], v[-1, 1], v[-1, 2]
        else:
            Qd = Qu = Ad = float("nan")
        Au = A - Ad
        if tU is not None and na > 0:
            Qd_l, Qu_l, Ad_l = lda_flow(RA, AN, Wmean[iz - 1], R)
        else:
            Qd_l = Qu_l = Ad_l = float("nan")
        rows.append((z, vi_m, vi_s, vm, Qd, Qu, Ad, Au,
                     Qd / Ad if Ad > 0 else float("nan"), Qu / Au if Au > 0 else float("nan"),
                     Qd_l, Qu_l, Qd_l / Ad_l if Ad_l > 0 else float("nan")))
    hdr = ("z,void_inst_mean,void_inst_std,void_meanField,Q_down,Q_up,A_down,A_up,u_down,u_up,"
           "Q_down_LDApts,Q_up_LDApts,u_down_LDApts")
    with open(os.path.join(out, "planes.csv"), "w") as f:
        f.write(hdr + "\n")
        for r in rows:
            f.write(",".join(f"{x:.6g}" for x in r) + "\n")

    # ---- monitors
    lines = [f"case: {case}", f"averaging window: t >= {tavg} s, probe samples: {nsamp}"]
    t, v = read_dat(os.path.join(pp, "holdup"), "volFieldValue.dat")
    if t is not None:
        s = t >= tavg if (t >= tavg).any() else t >= t[-1]
        lines.append(f"holdup (liquid region): mean {v[s,0].mean():.5f}  std {v[s,0].std():.5f}  last {v[-1,0]:.5f} at t={t[-1]:g}")
    t, v = read_dat(os.path.join(pp, "gasOut"), "surfaceFieldValue.dat")
    if t is not None:
        s = t >= tavg if (t >= tavg).any() else t >= t[-1]
        lines.append(f"gas out through top: mean {v[s,0]:.4e} kg/s" if False else
                     f"gas out through top: mean {v[s,0].mean():.4e} kg/s (last {v[-1,0]:.4e})")
    ls = time_dirs(os.path.join(pp, "freeSurfaceLines"))
    if ls:
        for name in ("axis", "rHalf"):
            p = glob.glob(os.path.join(ls[-1], f"{name}_alpha.water.xy"))
            if p:
                zz, aw = np.loadtxt(p[0], unpack=True)
                k = np.where((aw[:-1] >= 0.5) & (aw[1:] < 0.5))[0]
                if len(k):
                    k = k[-1]
                    zl = zz[k] + (0.5 - aw[k]) * (zz[k + 1] - zz[k]) / (aw[k + 1] - aw[k])
                    lines.append(f"free surface ({name}) at t={os.path.basename(ls[-1])}: z = {zl*1e3:.2f} mm, "
                                 f"swell = {(zl-cfg.LIQUID_HEIGHT)*1e3:.2f} mm")
    with open(os.path.join(out, "monitors.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")

    print("\n".join(lines))
    print("\n" + hdr)
    for r in rows:
        print(",".join(f"{x:.4g}" for x in r))
    print(f"\nwritten to {out}")

    if "--plot" in sys.argv:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(13, 4))
        for iz, z in enumerate(Z):
            ax[0].plot(np.array(RA) * 1e3, np.nanmean(Wmean[iz], 0), "o-", label=f"z = {z} m")
            ax[1].plot(np.array(RA) * 1e3, np.nanmean(Wstd[iz], 0), "o-")
            ax[2].plot(np.array(RA) * 1e3, np.nanmean(Amean[iz], 0), "o-")
        for a, yl in zip(ax, ("mean U_z liquid [m/s]", "sigma_z [m/s]", "mean alpha_air [-]")):
            a.set_xlabel("r [mm]"); a.set_ylabel(yl); a.grid(alpha=0.3)
        ax[0].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(out, "lda_radial.png"), dpi=150)
        print("plot:", os.path.join(out, "lda_radial.png"))


if __name__ == "__main__":
    main()
