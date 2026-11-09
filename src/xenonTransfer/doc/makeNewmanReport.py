#!/usr/bin/env python3
"""Build NewmanLiquidSideModel.pdf: liquid-side xenon transfer from a falling
molten-salt droplet (problem, analytical background, implementation plan,
verification plan, residence time).

Every number in the report is computed here. Run: python3 makeNewmanReport.py
Needs numpy, matplotlib, weasyprint.
"""

import math
import pathlib
import tempfile

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from weasyprint import HTML

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE/"NewmanLiquidSideModel.pdf"

# ---------------------------------------------------------------- inputs
d0 = 3.0e-4            # droplet diameter [m]
D_l = 2.006e-9         # Xe diffusivity in salt [m2/s] (thesis value)
rho_l = 4124.77        # salt density [kg/m3]
mu_l = 6.068e-3        # salt viscosity [Pa s]
T = 920.0              # [K]
p = 5000.0             # [Pa]
M_He = 4.0026e-3       # [kg/mol]
R_u = 8.314462618
rho_g = p*M_He/(R_u*T)  # helium density [kg/m3]
mu_g = 4.6913e-5       # helium viscosity [Pa s]
k_H = 2.08e-4          # dimensionless Henry (c_l/c_g at equilibrium)
D_g_thesis = 1.0e-5    # value used in the thesis case [m2/s]
g = 9.81
Y0 = 4.0e-5            # initial Xe mass fraction in the drop

# Measured in this work (xenonTransferInterFoam, run/xenonDroplet, 8 ranks, 0-0.3 ms)
MEAS = dict(
    areaRatio_t0=0.712, dp_t0=0.418e-3, beta_t0=0.138,
    areaRatio_late=1.14, dp_late=0.26e-3, areaUsed_late=0.61, beta_late=0.213,
    massBalance=0.9995, conservation=2.5e-10,
    # same run with the fixes of commit 12128de (perDropVol, volumeSphere, correctDropArea, shell slip)
    beta_fixed=0.2709, beta_fixed_spread=0.003, beta_nominal_fixed=0.272,
    areaScale_lo=1.40, areaScale_hi=1.60, massBalance_fixed=0.99996,
    Urel_iface=1.65e-3, Uslip_shell=7.9e-3, Ufall=1.56e-3, Umax=0.36,
)
THESIS = dict(t_end=0.0115, removed_Sh2=0.22, beta_fit_Sh2=0.191)

# ---------------------------------------------------------------- Newman solution
PI = math.pi
FO_SWITCH = 0.05


def newmanF(Fo):
    """Fraction released from a sphere, zero surface concentration (Crank eq. 6.20/6.21)."""
    Fo = np.atleast_1d(np.asarray(Fo, dtype=float))
    F = np.empty_like(Fo)
    s = Fo < FO_SWITCH
    F[s] = 6.0*np.sqrt(Fo[s]/PI) - 3.0*Fo[s]
    for i in np.where(~s)[0]:
        n = np.arange(1, 400)
        F[i] = 1.0 - 6.0/PI**2*np.sum(np.exp(-n**2*PI**2*Fo[i])/n**2)
    return F


def newmanSh(Fo):
    """Bulk-referred Sherwood number Sh = k_L d/D with k_L = (R/3)(dF/dt)/(1-F)."""
    Fo = np.atleast_1d(np.asarray(Fo, dtype=float))
    Sh = np.empty_like(Fo)
    s = Fo < FO_SWITCH
    x = Fo[s]
    Sh[s] = (2.0/3.0)*(3.0/np.sqrt(PI*x) - 3.0)/(1.0 - 6.0*np.sqrt(x/PI) + 3.0*x)
    for i in np.where(~s)[0]:
        n = np.arange(1, 400)
        e = np.exp(-n**2*PI**2*Fo[i])
        Sh[i] = (2.0*PI**2/3.0)*np.sum(e)/np.sum(e/n**2)
    return Sh


def newmanSeriesF(Fo, nmax=200000):
    """Series form only, used to check the short-time branch."""
    n = np.arange(1, nmax)
    return 1.0 - 6.0/PI**2*np.sum(np.exp(-n**2*PI**2*Fo)/n**2)


def newmanSeriesSh(Fo, nmax=200000):
    n = np.arange(1, nmax)
    e = np.exp(-n**2*PI**2*Fo)
    return (2.0*PI**2/3.0)*np.sum(e)/np.sum(e/n**2)


def FoForF(F):
    lo, hi = 1e-12, 5.0
    for _ in range(200):
        mid = math.sqrt(lo*hi)
        if newmanF(mid)[0] < F:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo*hi)


def wellMixedF(t, Sh, d=d0):
    return 1.0 - np.exp(-6.0*Sh*D_l*np.asarray(t)/d**2)


# Short-time branch vs. full series at the switch point (same Fo)
x = FO_SWITCH
switchErrF = abs((6*math.sqrt(x/PI) - 3*x) - newmanSeriesF(x))/newmanSeriesF(x)
switchErrSh = abs((2/3)*(3/math.sqrt(PI*x) - 3)/(1 - 6*math.sqrt(x/PI) + 3*x)
                  - newmanSeriesSh(x))/newmanSeriesSh(x)

# ---------------------------------------------------------------- gas side
def fullerDg(T, p_bar):
    MAB = 2.0/(1.0/4.0026 + 1.0/131.29)
    return 1e-4*0.00143*T**1.75/(p_bar*math.sqrt(MAB)*(2.67**(1/3) + 32.7**(1/3))**2)


def chapmanEnskogDg(T, p_bar):
    MAB = 2.0/(1.0/4.0026 + 1.0/131.29)
    sAB = 0.5*(2.551 + 4.047)
    Ts = T/math.sqrt(10.22*231.0)
    Om = (1.06036/Ts**0.15610 + 0.19300/math.exp(0.47635*Ts)
          + 1.03587/math.exp(1.52996*Ts) + 1.76474/math.exp(3.89411*Ts))
    return 1e-4*0.00266*T**1.5/(p_bar*math.sqrt(MAB)*sAB**2*Om)


D_g_F = fullerDg(T, p/1e5)
D_g_CE = chapmanEnskogDg(T, p/1e5)
D_g = 0.5*(D_g_F + D_g_CE)

kL_inf = (2*PI**2/3)*D_l/d0
kg_min = 2.0*D_g/d0
resRatio = (k_H/kg_min)/(1.0/kL_inf)
resRatio_thesisDg = (k_H/(2.0*D_g_thesis/d0))/(1.0/kL_inf)
eqRatio = k_H*rho_g/rho_l

# ---------------------------------------------------------------- free fall
lam = mu_g/p*math.sqrt(PI*R_u*T/(2*M_He))   # helium mean free path [m]


def cunningham(d):
    Kn = 2*lam/d
    return 1 + Kn*(1.257 + 0.4*math.exp(-1.1/Kn))


def fall(d, tEnd=20.0, slip=True):
    tau = rho_l*d**2/(18*mu_g)*(cunningham(d) if slip else 1.0)
    gEff = g*(1 - rho_g/rho_l)
    n = 200000
    t = np.geomspace(1e-6, tEnd, n)
    t = np.concatenate(([0.0], t))
    v = np.zeros_like(t)
    z = np.zeros_like(t)
    for i in range(1, len(t)):
        dt = t[i] - t[i-1]
        Re = rho_g*v[i-1]*d/mu_g
        f = 1 + 0.15*Re**0.687
        # implicit in drag, explicit in Re
        v[i] = (v[i-1] + gEff*dt)/(1 + f*dt/tau)
        z[i] = z[i-1] + 0.5*(v[i] + v[i-1])*dt
    return t, z, v, tau


def heightAt(tq, t, z):
    return float(np.interp(tq, t, z))


diams = [1.0e-4, 3.0e-4, 1.0e-3]
fallData = {d: fall(d) for d in diams}
tau0 = fallData[d0][3]
vterm0 = fallData[d0][2][-1]
Re_term0 = rho_g*vterm0*d0/mu_g
t5mm = float(np.interp(5e-3, fallData[d0][1], fallData[d0][0]))
z_thesis = heightAt(THESIS["t_end"], fallData[d0][0], fallData[d0][1])

targets = [0.5, 0.9, 0.99]
reqRows = []
for d in diams:
    R = d/2
    for Ft in targets:
        Fo = FoForF(Ft)
        t = Fo*R**2/D_l
        tt, zz, vv, _ = fallData[d]
        H = heightAt(t, tt, zz) if t <= tt[-1] else float("nan")
        reqRows.append((d, Ft, Fo, t, H))

# thesis-style extrapolation (first order with beta_fit, H = g t^2 / 2)
tF_thesis = {Ft: -math.log(1 - Ft)/THESIS["beta_fit_Sh2"] for Ft in targets}

# Internal circulation estimate (Hadamard-Rybczynski surface velocity)
kappa = mu_l/mu_g


def PeInt(U, d=d0):
    return U/(2*(1 + kappa))*d/D_l


def nCirc(z, d=d0):
    """Internal turnovers after a fall z: surface speed U/(2(1+kappa)) over a path pi*d/2."""
    return z/((1 + kappa)*PI*d)


U_thesis = float(np.interp(THESIS["t_end"], fallData[d0][0], fallData[d0][2]))
z50 = [r[4] for r in reqRows if r[0] == d0 and r[1] == 0.5][0]
z90 = [r[4] for r in reqRows if r[0] == d0 and r[1] == 0.9][0]


# ---------------------------------------------------------------- tables at thesis times
Fo_th = D_l*THESIS["t_end"]/(d0/2)**2
times = [1e-3, 5e-3, THESIS["t_end"], t5mm, 0.1, 1.0]
timeRows = []
for t in times:
    Fo = D_l*t/(d0/2)**2
    timeRows.append((t, Fo, newmanF(Fo)[0], wellMixedF(t, 2.0), wellMixedF(t, 2*PI**2/3),
                     newmanSh(Fo)[0], math.sqrt(PI*D_l*t)))

penRows = []
for Fo in [1e-4, 1e-3, 1e-2, 3e-2, 1e-1]:
    sh = newmanSh(Fo)[0]
    pen = 2/math.sqrt(PI*Fo)
    penRows.append((Fo, sh, pen, sh/pen))

# ---------------------------------------------------------------- figures
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e6e5e1"
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"   # reference palette slots 1-3 (light)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": INK2,
    "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2,
    "svg.fonttype": "path",
})

tmp = pathlib.Path(tempfile.mkdtemp())


def save(fig, name):
    path = tmp/f"{name}.svg"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path.as_uri()


# Fig 1: removal vs contact time
tt = np.geomspace(1e-4, 10.0, 400)
Fo = D_l*tt/(d0/2)**2
fig, ax = plt.subplots(figsize=(6.3, 3.4))
ax.plot(tt, 100*newmanF(Fo), color=C1)
ax.plot(tt, 100*wellMixedF(tt, 2*PI**2/3), color=C3)
ax.plot(tt, 100*wellMixedF(tt, 2.0), color=C2)
ax.plot([THESIS["t_end"]], [THESIS["removed_Sh2"]], "o", ms=7, mfc="white", mec=INK, mew=1.5)
ax.annotate("thesis CFD, Sh = 2\n(0.22 % at 11.5 ms)", (THESIS["t_end"], THESIS["removed_Sh2"]),
            xytext=(2.2e-4, 30), color=INK2, fontsize=8,
            arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
ax.text(3e-2, 45, "Newman, stagnant drop (exact)", color=INK, fontsize=8, ha="right")
ax.text(1.2, 100*wellMixedF(1.2, 2*PI**2/3) - 14, "well mixed, Sh = 6.58", color=INK, fontsize=8)
ax.text(1.2, 100*wellMixedF(1.2, 2.0) + 3, "well mixed, Sh = 2", color=INK, fontsize=8)
ax.axvline(t5mm, color=INK2, lw=0.8, ls=(0, (3, 3)))
ax.text(t5mm*1.08, 92, "5 mm fall", color=INK2, fontsize=8)
ax.set_xscale("log")
ax.set_xlim(1e-4, 10)
ax.set_ylim(0, 100)
ax.set_xlabel("contact time t (s)")
ax.set_ylabel("xenon removed from drop (%)")
fig1 = save(fig, "removal")

# Fig 2: Sh(Fo)
FoP = np.geomspace(1e-6, 1.0, 300)
fig, ax = plt.subplots(figsize=(6.3, 3.2))
ax.plot(FoP, newmanSh(FoP), color=C1)
ax.plot(FoP, 2/np.sqrt(PI*FoP), color=C2, ls=(0, (4, 2)))
ax.axhline(2*PI**2/3, color=C3, ls=(0, (4, 2)))
ax.text(2e-6, 2/math.sqrt(PI*2e-6)*1.25, "penetration 2/(π Fo)½", color=INK, fontsize=8)
ax.text(1.5e-2, 2*PI**2/3*0.62, "2π²/3 = 6.58", color=INK, fontsize=8)
ax.text(4e-2, newmanSh(4e-2)[0]*2.6, "Newman, bulk-referred", color=INK, fontsize=8, ha="center")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(1e-6, 1)
ax.set_xlabel(r"Fourier number $Fo = D_\ell\, t\, /\, R^2$")
ax.set_ylabel(r"$Sh_\ell = k_\ell\, d\, /\, D_\ell$")
ax.text(Fo_th, newmanSh(Fo_th)[0]*1.6, "11.5 ms", color=INK2, fontsize=8, ha="center")
ax.plot([Fo_th], newmanSh(Fo_th), "o", ms=6, mfc="white", mec=INK, mew=1.2)
fig2 = save(fig, "sherwood")

# Fig 3: removal vs fall height for three diameters
fig, ax = plt.subplots(figsize=(6.3, 3.4))
for d, c, lab, zl in zip(diams, [C3, C1, C2], ["d = 0.1 mm", "d = 0.3 mm", "d = 1 mm"], [0.04, 0.6, 0.6]):
    t, z, v, _ = fallData[d]
    m = z > 1e-5
    F = newmanF(D_l*t[m]/(d/2)**2)
    ax.plot(z[m], 100*F, color=c)
    Fl = 100*float(np.interp(zl, z[m], F))
    ax.text(zl, Fl + 7, lab, color=INK, fontsize=8, ha="right")
ax.axvline(5e-3, color=INK2, lw=0.8, ls=(0, (3, 3)))
ax.text(5.6e-3, 92, "thesis column (5 mm)", color=INK2, fontsize=8)
ax.set_xscale("log")
ax.set_xlim(1e-4, 10)
ax.set_ylim(0, 100)
ax.set_xlabel("free-fall height in helium H (m)")
ax.set_ylabel("xenon removed, Newman (%)")
fig3 = save(fig, "height")

# ---------------------------------------------------------------- equations (mathtext -> svg)
EQ = {}


def eq(name, tex, size=13):
    fig = plt.figure(figsize=(0.01, 0.01))
    fig.text(0, 0, f"${tex}$", fontsize=size, color=INK)
    path = tmp/f"eq_{name}.svg"
    fig.savefig(path, bbox_inches="tight", pad_inches=0.04, transparent=True)
    plt.close(fig)
    EQ[name] = f'<img class="eq" src="{path.as_uri()}"/>'


eq("model", r"\dot m''' = C\,K\,a_i\,\rho_\ell\,(Y_\ell - Y_{eq}),\qquad K=\left(\frac{1}{k_\ell}+\frac{k_H}{k_g}\right)^{-1},\qquad k_\ell=\frac{Sh_\ell\,D_\ell}{d_p}")
eq("yeq", r"Y_{eq}=k_H\,\frac{\rho_g}{\rho_\ell}\,Y_g")
eq("pde", r"\frac{\partial c}{\partial t}=\frac{D_\ell}{r^2}\frac{\partial}{\partial r}\left(r^2\frac{\partial c}{\partial r}\right),\quad c(r,0)=c_0,\quad c(R,t)=c_s\approx 0,\quad \left.\frac{\partial c}{\partial r}\right|_{0}=0")
eq("F", r"F(t)=1-\frac{\bar c(t)}{c_0}=1-\frac{6}{\pi^2}\sum_{n=1}^{\infty}\frac{1}{n^2}\exp\left(-n^2\pi^2 Fo\right),\qquad Fo=\frac{D_\ell\,t}{R^2}")
eq("Fshort", r"F\simeq 6\sqrt{\frac{Fo}{\pi}}-3\,Fo\qquad (Fo\ll 1)")
eq("kL", r"\frac{d\bar c}{dt}=-\frac{3}{R}\,k_\ell(t)\,(\bar c-c_s)\quad\Rightarrow\quad k_\ell(t)=\frac{R}{3}\,\frac{dF/dt}{1-F}")
eq("Sh", r"Sh_\ell(Fo)=\frac{k_\ell d}{D_\ell}=\frac{2\pi^2}{3}\,\frac{\sum_n e^{-n^2\pi^2Fo}}{\sum_n e^{-n^2\pi^2Fo}/n^2}\ \to\ \frac{2\pi^2}{3}=6.58\ \ (Fo\to\infty)")
eq("Shpen", r"Sh_\ell\simeq\frac{2}{\sqrt{\pi\,Fo}}\ \ \Leftrightarrow\ \ k_\ell=\sqrt{\frac{D_\ell}{\pi t}}\qquad(\mathrm{Higbie\ penetration},\ Fo\ll 1)")
eq("delta", r"\delta(t)\approx\sqrt{\pi D_\ell t}")
eq("kb", r"Sh_\infty^{KB}\approx 17.9,\qquad Pe_{int}=\frac{U\,d}{2(1+\kappa)\,D_\ell},\quad \kappa=\mu_\ell/\mu_g")
eq("fall", r"\frac{dv}{dt}=g\left(1-\frac{\rho_g}{\rho_\ell}\right)-\frac{v\,(1+0.15\,Re^{0.687})}{\tau_p},\qquad \tau_p=\frac{\rho_\ell d^2}{18\,\mu_g}\,C_c")
eq("age", r"\frac{\partial(\alpha_\ell\,\theta)}{\partial t}+\nabla\cdot(\alpha_\ell\,\phi\,\theta)=\alpha_\ell\,\chi_{drop}")
eq("renew", r"\theta_j^{n+1}=\left(\theta_j^{n}+\Delta t\right)\,\min\left(1,\frac{A_j^{n}}{A_j^{n+1}}\right)")
eq("flux", r"\dot M_j=k_\ell(Fo_j)\,A_j\,\rho_\ell\,\left(\bar Y_{\ell,j}-\bar Y_{eq,j}\right),\qquad Fo_j=\frac{4D_\ell\,\theta_j}{d_j^2}")


# ---------------------------------------------------------------- html helpers
def f(x, n=3):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "&ndash;"
    if x == 0:
        return "0"
    a = abs(x)
    if a >= 1e4 or a < 1e-3:
        m, e = f"{x:.{n-1}e}".split("e")
        return f"{m}&times;10<sup>{int(e)}</sup>"
    return f"{x:.{n}g}"


def table(head, rows, cls=""):
    h = "".join(f"<th>{c}</th>" for c in head)
    b = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f'<table class="{cls}"><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table>'


paramRows = [
    ("Drop diameter d", "0.30 mm", "thesis"),
    ("Xe diffusivity in salt D<sub>&#8467;</sub>", f"{f(D_l)} m&sup2;/s", "thesis value, source to be checked"),
    ("Salt density / viscosity", f"{rho_l:.0f} kg/m&sup3; / {mu_l*1e3:.2f} mPa&middot;s", "thesis"),
    ("He density at 920 K, 5 kPa", f"{f(rho_g)} kg/m&sup3;", "ideal gas"),
    ("He viscosity", f"{f(mu_g)} Pa&middot;s", "thesis"),
    ("He&ndash;Xe diffusivity D<sub>g</sub>", f"{f(D_g_F)} (Fuller) / {f(D_g_CE)} (Chapman&ndash;Enskog) m&sup2;/s", "[13, 14, 15]"),
    ("Dimensionless Henry coefficient k<sub>H</sub>", f"{f(k_H)}", "thesis; cf. [16, 17]"),
    ("He mean free path / Knudsen number 2&lambda;/d", f"{f(lam*1e6)} &micro;m / {f(2*lam/d0)}", "kinetic theory"),
]

resistRows = [
    ("Liquid side 1/k<sub>&#8467;</sub> (Sh = 6.58)", f"{f(1/kL_inf)} s/m"),
    ("Gas side k<sub>H</sub>/k<sub>g</sub> (Sh<sub>g</sub> = 2, physical D<sub>g</sub>)", f"{f(k_H/kg_min)} s/m"),
    ("ratio gas/liquid", f"{f(resRatio)}"),
    ("ratio with the thesis D<sub>g</sub> = 10<sup>&minus;5</sup>", f"{f(resRatio_thesisDg)}"),
    ("Equilibrium back-pressure Y<sub>eq</sub>/Y<sub>g</sub> = k<sub>H</sub>&rho;<sub>g</sub>/&rho;<sub>&#8467;</sub>", f"{f(eqRatio)}"),
]

timeTab = table(
    ["t", "Fo", "Newman F", "Sh = 2 (thesis)", "Sh = 6.58", "Sh<sub>&#8467;</sub>(Fo)", "&delta; = (&pi;D<sub>&#8467;</sub>t)<sup>&frac12;</sup>"],
    [(f"{t*1e3:.4g} ms", f(Fo), f"{100*F:.3g} %", f"{100*a:.3g} %", f"{100*b:.3g} %", f"{sh:.3g}", f"{dl*1e6:.3g} &micro;m")
     for t, Fo, F, a, b, sh, dl in timeRows])

penTab = table(["Fo", "Sh<sub>&#8467;</sub> Newman", "2/(&pi;Fo)<sup>&frac12;</sup>", "ratio"],
               [(f(Fo), f"{a:.4g}", f"{b:.4g}", f"{c:.3f}") for Fo, a, b, c in penRows])

reqTab = table(
    ["d", "target F", "Fo", "contact time", "free-fall height", "thesis extrapolation (Sh = 2, &frac12;gt&sup2;)"],
    [(f"{d*1e3:.2g} mm", f"{100*Ft:.0f} %", f(Fo), f"{t:.3g} s", (f"{H:.3g} m" if not math.isnan(H) else "&gt; 20 s of fall"),
      (f"{tF_thesis[Ft]:.3g} s / {f(0.5*g*tF_thesis[Ft]**2)} m" if d == d0 else ""))
     for d, Ft, Fo, t, H in reqRows])

fallTab = table(
    ["d", "&tau;<sub>p</sub>", "terminal velocity", "Re at terminal velocity", "t to fall 5 mm"],
    [(f"{d*1e3:.2g} mm", f"{fallData[d][3]:.3g} s", f"{fallData[d][2][-1]:.3g} m/s",
      f"{rho_g*fallData[d][2][-1]*d/mu_g:.3g}", f"{float(np.interp(5e-3, fallData[d][1], fallData[d][0]))*1e3:.3g} ms")
     for d in diams])

measRows = [
    ("Iso-surface area / sphere area", f"{MEAS['areaRatio_t0']:.2f}", f"{MEAS['areaRatio_late']:.2f}"),
    ("d<sub>p</sub> = 6V/A (true 0.30 mm)", f"{MEAS['dp_t0']*1e3:.3g} mm", f"{MEAS['dp_late']*1e3:.2g} mm"),
    ("Fraction of iso area that carries flux", "1.00", f"{MEAS['areaUsed_late']:.2f}"),
    ("Drop removal rate &beta; (nominal 6ShD/d&sup2; = 0.267 s<sup>&minus;1</sup>)", f"{MEAS['beta_t0']:.3g} s<sup>&minus;1</sup>", f"{MEAS['beta_late']:.3g} s<sup>&minus;1</sup>"),
]

fixedNote = (
    f"<p>With the fixes (<code>dpMode perDropVol</code>, <code>dropArea volumeSphere</code>, <code>correctDropArea true</code>, "
    f"<code>slipVelocity shell</code>) the same run gives a drop removal rate, measured from the liquid xenon inventory, of "
    f"&beta; = {MEAS['beta_fixed']:.4g} s<sup>&minus;1</sup> against the nominal {MEAS['beta_nominal_fixed']} s<sup>&minus;1</sup> "
    f"(6&middot;Sh&middot;D<sub>&#8467;</sub>/D<sub>vol</sub>&sup2;), constant to {100*MEAS['beta_fixed_spread']:.1g} % while the interface "
    f"smears; the flux cells are rescaled by {MEAS['areaScale_lo']:.2f}&ndash;{MEAS['areaScale_hi']:.2f}. Liquid loss / integrated model "
    f"rate = {MEAS['massBalance_fixed']}.</p>"
    f"<div class=\"warn\"><b>Parasitic currents.</b> At 0.3 ms the maximum velocity in the domain is {MEAS['Umax']} m/s, in the helium "
    f"just above the drop, while the drop falls at {f(MEAS['Ufall'])} m/s. The capillary velocity &sigma;/&mu;<sub>g</sub> is about "
    f"{0.2617/mu_g:.0f} m/s and the density ratio about {f(rho_l/rho_g, 2)}, so spurious currents in the gas are expected. As a consequence the "
    f"shell slip velocity is {f(MEAS['Uslip_shell'])} m/s (the old adjacent-cell value was {f(MEAS['Urel_iface'])} m/s): the shell "
    f"definition is correct, but its input is polluted at these low drop velocities. This does not affect the present results (the slip only enters "
    f"the negligible gas side) but it must be reduced before any velocity-dependent closure (internal circulation, Higbie, Kronig&ndash;Brink "
    f"switching) is trusted: a static-drop test (g = 0) quantifies it directly.</div>")

# ---------------------------------------------------------------- document
css = """
@page { size: A4; margin: 20mm 19mm 20mm 19mm;
        @bottom-center { content: counter(page); font-size: 8.5pt; color: #52514e; } }
body { font-family: 'DejaVu Serif', serif; font-size: 10pt; line-height: 1.42; color: #0b0b0b; }
h1 { font-family: 'DejaVu Sans'; font-size: 17pt; margin: 0 0 2mm 0; }
h2 { font-family: 'DejaVu Sans'; font-size: 12.5pt; margin: 7mm 0 2mm 0; border-bottom: 0.6pt solid #c9c8c2; padding-bottom: 1mm; }
h3 { font-family: 'DejaVu Sans'; font-size: 10.5pt; margin: 4.5mm 0 1.5mm 0; }
.sub { color: #52514e; font-family: 'DejaVu Sans'; font-size: 9.5pt; margin-bottom: 5mm; }
.eq { display: block; margin: 2.5mm auto; max-width: 100%; }
figure { margin: 3mm 0 4mm 0; text-align: center; page-break-inside: avoid; }
figure img { width: 100%; }
figcaption { font-size: 8.8pt; color: #52514e; text-align: left; margin-top: 1mm; }
table { border-collapse: collapse; width: 100%; font-size: 8.6pt; margin: 2.5mm 0 3.5mm 0; page-break-inside: avoid; }
th { font-family: 'DejaVu Sans'; text-align: left; font-weight: bold; border-bottom: 0.8pt solid #0b0b0b; padding: 1.2mm 1.6mm; }
td { border-bottom: 0.4pt solid #e1e0db; padding: 1mm 1.6mm; vertical-align: top; }
.box { border-left: 2.2pt solid #2a78d6; background: #f4f7fc; padding: 2.5mm 4mm; margin: 3mm 0; }
.warn { border-left: 2.2pt solid #eb6834; background: #fdf5f1; padding: 2.5mm 4mm; margin: 3mm 0; }
code { font-family: 'DejaVu Sans Mono'; font-size: 8.6pt; }
ol.refs li { font-size: 8.8pt; margin-bottom: 1.2mm; }
"""

html = f"""<!DOCTYPE html><html><head><meta charset="utf-8"><style>{css}</style></head><body>

<h1>Liquid-side xenon transfer from a falling molten-salt droplet</h1>
<div class="sub">Problem statement, analytical background (Newman), implementation plan, verification plan and residence-time
analysis for <code>xenonTransferInterFoam</code> / <code>gasDiffusionSimple</code> (tom-v2606). All numbers in this document are
computed by <code>src/xenonTransfer/doc/makeNewmanReport.py</code>.</div>

<div class="box"><b>Summary.</b> In the falling-droplet case the gas-side resistance and the equilibrium back-pressure are
negligible (ratios of order 10<sup>&minus;10</sup>), so the problem reduces to diffusion of xenon out of a sphere whose surface concentration
is zero: Newman&rsquo;s problem [1, 2]. The thesis closes the liquid side with a constant Sh<sub>&#8467;</sub> = 2, which is the
<i>external</i> steady value for a sphere in a stagnant medium and has no meaning for transport <i>inside</i> a drop. The exact
inner solution removes {100*timeRows[2][2]:.1f} % of the xenon in the 11.5 ms simulated in the thesis, against
{THESIS['removed_Sh2']} % obtained with Sh<sub>&#8467;</sub> = 2: about {timeRows[2][2]*100/THESIS['removed_Sh2']:.0f}&times; more. The concentration
boundary layer inside the drop is a few micrometres thick at these times and cannot be resolved by the VOF mesh (16&ndash;30 &micro;m cells),
so a subgrid closure k<sub>&#8467;</sub>(t) is required. This document derives it, proposes how to carry the contact time
through droplet series, re-entry and splashing, and lists what must be verified before using it.</div>

<h2>1. Problem</h2>
<h3>1.1 Current model</h3>
<p><code>gasDiffusionSimple</code> computes a volumetric xenon sink in the interface cells of each identified drop:</p>
{EQ['model']}
{EQ['yeq']}
<p>with a<sub>i</sub> the interface area density from the &alpha; = 0.5 iso-surface, Y<sub>&#8467;</sub> the liquid xenon mass
fraction sampled next to the interface, and Y<sub>g</sub> the gas-side value across the interface. The liquid side is closed with a
<code>fixed</code> Sh<sub>&#8467;</sub> (2 in the thesis reference case, 2&ndash;100 in the sweep) or with the Higbie form used for bubbles.</p>

<h3>1.2 Gas side and back-pressure are negligible</h3>
{table(["quantity", "value"], resistRows)}
<p style="font-size:8.6pt;color:#52514e">Table 1. Resistances in series for d = 0.3 mm (series resistances 1/k<sub>&#8467;</sub> and
k<sub>H</sub>/k<sub>g</sub> of K) and equilibrium back-pressure. Inputs:</p>
{table(["parameter", "value", "source"], paramRows)}
<p>The gas side therefore drops out of K for any plausible Sh<sub>g</sub>, and the liquid sees a surface concentration of essentially
zero. Note that the thesis used D<sub>g</sub> = 10<sup>&minus;5</sup> m&sup2;/s; the physical He&ndash;Xe value at 920 K and 5 kPa is
about {f(D_g)} m&sup2;/s ({f(D_g_F)} from Fuller&ndash;Schettler&ndash;Giddings [13], {f(D_g_CE)} from Chapman&ndash;Enskog with
Lennard-Jones parameters and the Neufeld collision integral [14, 15]). This changes nothing in K, but the same value enters the resolved
gas-phase species diffusion, which is explicit in the code: a physical D<sub>g</sub> would require &Delta;t &lt; &Delta;x&sup2;/(6D<sub>g</sub>)
&asymp; 10<sup>&minus;8</sup> s and is not protected by <code>maxDi</code> (multiphaseInter returns a zero diffusion number).</p>

<h3>1.3 What Sh<sub>&#8467;</sub> = 2 means, and why it is wrong here</h3>
<p>Sh = 2 is the steady solution for diffusion from a sphere into an <i>unbounded</i> stagnant medium (the outside problem). For
transport <i>inside</i> a sphere the steady value does not exist (the drop is depleted), and the long-time limit of the
bulk-referred coefficient is 2&pi;&sup2;/3 = 6.58 [1, 2, 4]; at short times the coefficient is much larger because the depleted layer is
thin (penetration regime, Higbie [3]). The values 6.58 and 17.9 in the thesis sweep are indeed the stagnant (Newman) and the
circulating (Kronig&ndash;Brink [5]) long-time limits, but the transient was not considered.</p>

<h3>1.4 Why the CFD cannot resolve it</h3>
<p>The depleted layer grows as {EQ['delta']}
which is {timeRows[0][6]*1e6:.2g} &micro;m after 1 ms, {timeRows[2][6]*1e6:.2g} &micro;m after 11.5 ms and {timeRows[4][6]*1e6:.2g} &micro;m after
0.1 s. The droplet mesh has 16&ndash;30 &micro;m cells; resolving the layer with a few cells needs sub-micrometre spacing, i.e. hundreds of cells per
diameter in 3D. In addition the xenon storage Z = &alpha;Y is advected with first-order upwind, whose numerical diffusivity
(~U&Delta;x/2 &asymp; 10<sup>&minus;7</sup>&ndash;10<sup>&minus;6</sup> m&sup2;/s while the drop moves) exceeds D<sub>&#8467;</sub> by two to three orders
of magnitude, so the drop interior is numerically well mixed. The thesis D<sub>&#8467;</sub> sweep confirms it: the fitted removal rate scales
exactly linearly with D<sub>&#8467;</sub> (0.191, 1.91, 18.9, 183 s<sup>&minus;1</sup>), i.e. D<sub>&#8467;</sub> only acts through the input k<sub>&#8467;</sub>.
Unresolved interfacial boundary layers in VOF are a known problem and are normally handled with subgrid-scale closures [10, 11, 12].
The consistent formulation is therefore: <b>well-mixed drop + subgrid k<sub>&#8467;</sub>(t) driven by the drop-bulk concentration</b>.</p>

<h2>2. Analytical background</h2>
<h3>2.1 Diffusion out of a sphere (Newman 1931)</h3>
<p>For constant D<sub>&#8467;</sub>, no internal motion, uniform initial concentration c<sub>0</sub> and constant surface concentration
c<sub>s</sub> (here &asymp; 0):</p>
{EQ['pde']}
<p>the fraction released is [1; 2, eq. 6.20]</p>
{EQ['F']}
<p>with the short-time form [2, eq. 6.21; for Fo &lt; {FO_SWITCH} it differs from the full series by less than {f(switchErrF, 2)} in F
and {f(switchErrSh, 2)} in the Sh below]</p>
{EQ['Fshort']}

<h3>2.2 Equivalent bulk-referred coefficient</h3>
<p>Writing the drop balance with a coefficient referred to the bulk (volume-averaged) concentration, which is what a well-mixed CFD drop
carries, gives an exact time-dependent k<sub>&#8467;</sub>:</p>
{EQ['kL']}
{EQ['Sh']}
<p>At short times this tends to the penetration result</p>
{EQ['Shpen']}
<p>Integrating the drop balance with this k<sub>&#8467;</sub>(t) reproduces F(t) exactly, which is the key property for verification. The
penetration form is within 2 % of Newman for Fo &le; 10<sup>&minus;4</sup> and within 6 % for Fo &le; 10<sup>&minus;3</sup> (11.5 ms for
d = 0.3 mm), but 19 % low at Fo = 0.01 and a factor two low at Fo = 0.1, where depletion of the sphere interior matters (table 2). This matters for
splashing and films (section 3.2), where the sphere geometry is lost but the surface-age concept survives only in the penetration regime.</p>
{penTab}
<p style="font-size:8.6pt;color:#52514e">Table 2. Newman bulk-referred Sh vs. the penetration limit.</p>

<figure><img src="{fig2}"/><figcaption>Figure 1. Bulk-referred liquid-side Sherwood number of a stagnant sphere. Marker: the end of the thesis
simulation (Fo = {f(Fo_th)}, Sh = {newmanSh(Fo_th)[0]:.3g}).</figcaption></figure>

<h3>2.3 Consequences for the thesis case</h3>
{timeTab}
<p style="font-size:8.6pt;color:#52514e">Table 3. Removed fraction for d = 0.3 mm. &ldquo;Sh = 2&rdquo; and &ldquo;Sh = 6.58&rdquo;: well-mixed drop
with constant coefficient, F = 1 &minus; exp(&minus;6 Sh D<sub>&#8467;</sub> t/d&sup2;). The thesis CFD obtained 0.22 % at 11.5 ms with Sh = 2
(slightly below 0.31 % because of the area/diameter reconstruction issue of section 5).</p>
<figure><img src="{fig1}"/><figcaption>Figure 2. Xenon removed from a 0.3 mm drop vs. contact time: exact stagnant-drop solution and well-mixed
models with constant Sh.</figcaption></figure>

<h3>2.4 Validity conditions (to be checked, not assumed)</h3>
<ul>
<li><b>No internal circulation.</b> For a drop in a gas the viscosity ratio is &kappa; = &mu;<sub>&#8467;</sub>/&mu;<sub>g</sub> = {kappa:.0f},
and the Hadamard&ndash;Rybczynski interface velocity of a clean drop is U/(2(1+&kappa;)) [6, 7]. The internal P&eacute;clet number is
{EQ['kb']}
and the number of internal turnovers after a fall z is N &asymp; z/((1+&kappa;)&pi;d). At the end of the thesis window (U = {U_thesis:.2g} m/s)
Pe<sub>int</sub> = {PeInt(U_thesis):.2g} but N = {f(nCirc(z_thesis), 2)}; after 5 mm N = {nCirc(5e-3):.2g}; after the {z50:.2g} m needed for
50 % removal N = {nCirc(z50):.2g}; after {z90:.2g} m (90 %) N = {nCirc(z90):.2g}. Circulation is therefore irrelevant for millimetre falls,
but for decimetre-to-metre contactors it can develop and push transfer towards the Kronig&ndash;Brink limit (Sh<sub>&infin;</sub> &asymp; 17.9,
effective diffusivity 2.25 D<sub>&#8467;</sub> [5, 7]). Newman is the lower bound and Kronig&ndash;Brink the laminar upper bound; surface-active
impurities push towards Newman. The solver already prints the resolved internal-circulation velocity and P&eacute;clet number
(<code>debugInternalCirculation</code>), but see the parasitic-current caveat in section 5.</li>
<li><b>Uniform initial concentration</b> and <b>constant D<sub>&#8467;</sub></b> (dilute xenon, isothermal drop).</li>
<li><b>Near-spherical drop.</b> For deformed drops the area enters through A<sub>j</sub>; in the penetration regime the flux per unit area
does not depend on the shape.</li>
<li><b>Zero surface concentration</b>, verified in section 1.2.</li>
<li><b>Rarefaction.</b> The helium mean free path is {lam*1e6:.2g} &micro;m (Knudsen number 2&lambda;/d = {2*lam/d0:.2f}); this only
affects the drag (Cunningham factor {cunningham(d0):.3f}), not the liquid side.</li>
</ul>

<h2>3. Implementation plan (not yet implemented)</h2>
<h3>3.1 Closure</h3>
<p>New option <code>liquidSherwoodModel Newman;</code> in <code>gasDiffusionSimple</code>. For each drop j:</p>
{EQ['flux']}
<p>where &theta;<sub>j</sub> is the contact time of drop j, A<sub>j</sub> its transfer area (already available:
<code>dropArea iso | volumeSphere</code> with <code>correctDropArea</code>), and Y&#772;<sub>&#8467;,j</sub> the mass-weighted bulk xenon
fraction of the drop (already computed as <code>dropYbulk</code>). The drop flux is distributed over the drop&rsquo;s flux cells
proportionally to their area, so the existing conservative face-based mapping to the gas is reused unchanged. Sh(Fo) is evaluated
with the short-time expression for Fo &lt; {FO_SWITCH} and the series (truncated where n&sup2;&pi;&sup2;Fo &gt; 50) above.</p>
<p>Driving the flux with the bulk value instead of the interface-cell value removes the mesh dependence of the sampled
concentration and makes the CFD drop reproduce the Newman curve exactly in the limit where it applies.</p>

<h3>3.2 Contact time for single drops, series, re-entry and splashing</h3>
<p>The target applications (a series of droplets, the same droplet re-entering the domain after leaving it, splashing on an inclined
plate) rule out a global clock and make Lagrangian drop-identity tracking fragile (drops merge, split, leave and re-enter, and identities
are not preserved by the region-split clustering between steps or across processors). The proposal is an Eulerian <b>contact-age</b> field
&theta; carried by the liquid, transported with the same conservative Z = &alpha;&theta; machinery as xenon (mass-consistent with &alpha;):</p>
{EQ['age']}
<p>with &chi;<sub>drop</sub> = 1 in cells belonging to an identified drop (or film fragment) in contact with helium and 0 in the feed/injector
liquid. This is the &ldquo;local mean age&rdquo; concept of residence-time distributions [8, 9]. The drop value is the mass-weighted mean
&theta;<sub>j</sub>. Behaviour in the target cases:</p>
<ul>
<li><b>Series of drops:</b> each drop starts with &theta; = 0 at the injector and accumulates its own age; no identity bookkeeping.</li>
<li><b>Re-entry:</b> if the outlet is mapped back to the inlet, both Z<sub>Xe</sub> and Z<sub>&theta;</sub> are carried, so the drop resumes
with its own depletion history (F continuous across the boundary).</li>
<li><b>Breakup:</b> fragments inherit the parent age (exact to leading order in the penetration regime, where the depleted layer stays on
the surface).</li>
<li><b>Coalescence:</b> mass-weighted mean age; this is an approximation because k<sub>&#8467;</sub>(&theta;) is non-linear (state it as such).</li>
<li><b>Splashing / film formation:</b> new interface exposes undepleted liquid (surface renewal, Danckwerts [4]). Proposed per-drop renewal
when the drop area grows: {EQ['renew']}
i.e. the newly created area is assigned zero age and mixed by area with the old surface. On the plate the sphere geometry is lost; the
penetration form k<sub>&#8467;</sub> = (D<sub>&#8467;</sub>/&pi;&theta;)<sup>&frac12;</sup> remains valid while the depleted depth
&delta; = (&pi;D<sub>&#8467;</sub>&theta;)<sup>&frac12;</sup> is small compared with the film thickness; beyond that a film model is needed.
This renewal rule is a modelling hypothesis and must be validated, not only verified.</li>
</ul>

<h3>3.3 Numerical aspects</h3>
<ul>
<li>The source stays explicit: &beta;&Delta;t = 6Sh D<sub>&#8467;</sub>&Delta;t/d&sup2; &lt; 10<sup>&minus;3</sup> even at Fo = 10<sup>&minus;6</sup>.</li>
<li>Per-drop quantities need global reductions; for droplet series the ~20 scalar reductions per drop per step should be merged into a
single list reduction.</li>
<li>Sh(Fo) is singular as Fo &rarr; 0; the first step after release should use the time-averaged coefficient over the step,
k&#772; = (F(&theta;+&Delta;t) &minus; F(&theta;))/(&Delta;t &middot; 3(1&minus;F(&theta;))/R), which is exact for the step.</li>
</ul>

<h2>4. Verification</h2>
<h3>4.1 Why</h3>
<p>The closure replaces physics the mesh cannot resolve. The CFD then only supplies geometry (area, volume), motion and bookkeeping
(mass conservation, mapping to the gas, age transport). Verification checks that this chain reproduces the analytical solution exactly where
that solution applies, so that later differences in the target cases (series, re-entry, splashing) can be attributed to physics rather than to
implementation or discretisation. It is distinct from validation against experiments, which concerns whether Newman/renewal is the right
physics (e.g. internal circulation, section 2.4).</p>
<h3>4.2 What to verify</h3>
{table(["#", "test", "checks", "acceptance"], [
    ("V0", "Unit test of Sh(Fo), F(Fo)", "branch switch, limits 2/(&pi;Fo)<sup>&frac12;</sup> and 6.58, against an independent series", "relative error &lt; 10<sup>&minus;6</sup>"),
    ("V1", "0D integration of the drop balance with the solver time stepping", "explicit source and first-step treatment reproduce F(t)", "error &lt; 0.1 % in F up to Fo = 0.1"),
    ("V2", "Stagnant drop in quiescent helium (g = 0), three meshes", "area, bulk sampling, mapping to the gas, conservation, mesh independence", "F within 2 % of Newman at 10, 50, 100 ms on all meshes; liquid loss = &int;&#7745; within 0.1 %"),
    ("V3", "Falling drop (this case) to t &asymp; 35 ms (5 mm)", "insensitivity to motion and slip; internal-circulation diagnostics", "same F as V2 within 2 %; Pe<sub>int</sub> reported"),
    ("V4", "Two drops released at different times", "age transport and per-drop &theta;<sub>j</sub>", "&theta;<sub>j</sub> equals each drop&rsquo;s own contact time within one step"),
    ("V5", "Drop leaving the bottom and re-entering at the top (mapped boundary)", "conservation of Z<sub>Xe</sub>, Z<sub>&theta;</sub> across the mapping", "F(t) continuous, total Xe conserved to 10<sup>&minus;8</sup>"),
    ("V6", "Prescribed stretching of a drop (kinematic test)", "area-renewal rule", "&theta;<sub>j</sub> follows the renewal rule exactly"),
])}
<p style="font-size:8.6pt;color:#52514e">Table 6. Verification tests, in the order they should be run.</p>
<p>V2 is the central test: with g = 0 the Newman solution is exact for the CFD set-up, so any discrepancy is an implementation or
geometry error. Before the Newman model, the same test with a fixed Sh already checks the geometry chain (section 5).</p>

<h2>5. Geometry and slip fixes already applied (for reference)</h2>
<p>Measured on <code>run/xenonDroplet</code> with the original code (Sh = 2):</p>
{table(["quantity", "t = 0 (sharp interface)", "t &ge; 40 &micro;s (smeared)"], measRows)}
<p>The flux scales as k<sub>&#8467;</sub>A &prop; A&sup2;/V, so area errors enter squared. The face-based mapping to the gas is conservative
(liquid loss / integrated model rate = {MEAS['massBalance']}, total xenon conserved to {f(MEAS['conservation'])}). Fixes now in the code:
<code>dpMode perDropVol</code> (diameter from the conserved volume), <code>dropArea volumeSphere|iso</code> with <code>correctDropArea</code>
(each drop transfers through its full area), <code>slipVelocity shell</code> (drop velocity minus helium velocity averaged in a shell 1.5&ndash;3 radii
away; the phases share one velocity field, so the old adjacent-cell difference was not a slip velocity), and Y = Z/&alpha; reconstructed with the
start-of-step &alpha;.</p>
{fixedNote}

<h2>6. Residence time and fall height</h2>
<h3>6.1 Drop kinematics</h3>
<p>The drop is released from rest. With Schiller&ndash;Naumann drag and the Cunningham slip correction C<sub>c</sub>:</p>
{EQ['fall']}
{fallTab}
<p style="font-size:8.6pt;color:#52514e">Table 4. Fall kinematics in helium at 920 K, 5 kPa.</p>
<p>Because &tau;<sub>p</sub> = {tau0:.2g} s for d = 0.3 mm, the drop is in practically free fall over the thesis window: it falls only
{z_thesis*1e3:.2g} mm in 11.5 ms, so the 5 mm column is never traversed; reaching 5 mm takes {t5mm*1e3:.3g} ms. In the thesis set-up the drop
also starts 0.35 mm (about 12 cells) below the top boundary. Since transfer is controlled by the liquid side and there is no internal circulation in
this window, what matters is the contact time, not the height. Note that the drop never reaches its terminal velocity in any practical column
(&tau;<sub>p</sub> = {tau0:.2g} s, terminal velocity {vterm0:.2g} m/s at Re = {Re_term0:.2g}), so the free-fall height grows almost as
&frac12;gt&sup2; up to a few tenths of a second.</p>
<h3>6.2 Required contact time and height</h3>
{reqTab}
<p style="font-size:8.6pt;color:#52514e">Table 5. Contact time from Newman and the corresponding free-fall height (drag and slip included). Last column:
the thesis back-up estimate, i.e. first-order decay with the fitted &beta; = {THESIS['beta_fit_Sh2']} s<sup>&minus;1</sup> (Sh = 2) and drag-free fall
H = &frac12;gt&sup2;.</p>
<figure><img src="{fig3}"/><figcaption>Figure 3. Xenon removed vs. free-fall height in helium for three drop diameters (Newman, stagnant drop).
Small drops strip faster (Fo &prop; t/d&sup2;) but also fall slower.</figcaption></figure>
<div class="warn">The thesis extrapolation underestimates removal at short contact times by more than an order of magnitude (wrong Sh and an
exponential instead of a &radic;t law) and assumes drag-free fall at all heights. Conversely, at long contact times possible internal
circulation would increase transfer above Newman. Contactor heights should be estimated with the Newman curve as a conservative
reference and bracketed with Kronig&ndash;Brink once the circulation diagnostics are available.</div>

<h2>References</h2>
<p style="font-size:8.6pt;color:#52514e">Bibliographic details were written from memory and should be checked against the originals before citing.</p>
<ol class="refs">
<li>A. B. Newman, The drying of porous solids: diffusion and surface emission equations, Trans. AIChE 27 (1931) 203&ndash;220.</li>
<li>J. Crank, The Mathematics of Diffusion, 2nd ed., Clarendon Press, Oxford, 1975 (ch. 6, sphere; eqs. 6.18&ndash;6.21).</li>
<li>R. Higbie, The rate of absorption of a pure gas into a still liquid during short periods of exposure, Trans. AIChE 31 (1935) 365&ndash;389.</li>
<li>P. V. Danckwerts, Significance of liquid-film coefficients in gas absorption, Ind. Eng. Chem. 43 (1951) 1460&ndash;1467.</li>
<li>R. Kronig, J. C. Brink, On the theory of extraction from falling droplets, Appl. Sci. Res. A2 (1950) 142&ndash;154.</li>
<li>A. E. Handlos, T. Baron, Mass and heat transfer from drops in liquid&ndash;liquid extraction, AIChE J. 3 (1957) 127&ndash;136.</li>
<li>R. Clift, J. R. Grace, M. E. Weber, Bubbles, Drops, and Particles, Academic Press, New York, 1978 (ch. 3 and 5).</li>
<li>P. V. Danckwerts, Continuous flow systems: distribution of residence times, Chem. Eng. Sci. 2 (1953) 1&ndash;13.</li>
<li>F. Ghirelli, B. Leckner, Transport equation for the local residence time of a fluid, Chem. Eng. Sci. 59 (2004) 513&ndash;523.</li>
<li>D. Bothe, S. Fleckenstein, A volume-of-fluid-based method for mass transfer processes at fluid particles, Chem. Eng. Sci. 101 (2013) 283&ndash;302.</li>
<li>D. Deising, H. Marschall, D. Bothe, A unified single-field model framework for volume-of-fluid simulations of interfacial species transfer applied to bubbly flows, Chem. Eng. Sci. 139 (2016) 173&ndash;195.</li>
<li>A. Weiner, D. Bothe, Advanced subgrid-scale modeling for convection-dominated species transport at fluid interfaces with application to mass transfer from rising bubbles, J. Comput. Phys. 347 (2017) 261&ndash;289.</li>
<li>E. N. Fuller, P. D. Schettler, J. C. Giddings, A new method for prediction of binary gas-phase diffusion coefficients, Ind. Eng. Chem. 58(5) (1966) 18&ndash;27.</li>
<li>B. E. Poling, J. M. Prausnitz, J. P. O&rsquo;Connell, The Properties of Gases and Liquids, 5th ed., McGraw-Hill, 2001 (ch. 11; Lennard-Jones parameters, app. B).</li>
<li>T. R. Marrero, E. A. Mason, Gaseous diffusion coefficients, J. Phys. Chem. Ref. Data 1 (1972) 3&ndash;118 (measured He&ndash;Xe data to check the estimates).</li>
<li>M. Blander, W. R. Grimes, N. V. Smith, G. M. Watson, Solubility of noble gases in molten fluoride mixtures, J. Phys. Chem. 63 (1959) 1164&ndash;1167.</li>
<li>G. M. Watson, R. B. Evans, W. R. Grimes, N. V. Smith, Solubility of noble gases in molten fluorides. In lithium-beryllium fluoride, J. Chem. Eng. Data 7 (1962) 285&ndash;287.</li>
<li>R. J. Kedl, A. Houtzeel, Development of a model for computing <sup>135</sup>Xe migration in the MSRE, ORNL-4069, 1967.</li>
<li>E. F. Tiso, Numerical analysis of ex-core xenon removal from molten salt reactor fuel with a multicomponent volume-of-fluid solver, MSc thesis, Politecnico di Milano, 2026.</li>
</ol>
</body></html>"""

HTML(string=html, base_url=str(tmp)).write_pdf(OUT)
print(f"wrote {OUT}")
print(f"D_g Fuller={D_g_F:.4g} CE={D_g_CE:.4g}  rho_g={rho_g:.5g}  tau_p={tau0:.3g}  t5mm={t5mm:.4g}  z(11.5ms)={z_thesis:.4g}")
print(f"switch errors F={switchErrF:.2e} Sh={switchErrSh:.2e}")
for r in timeRows:
    print("t=%.4g Fo=%.3g Newman=%.4g Sh2=%.4g Sh658=%.4g Sh(Fo)=%.4g delta=%.3g" % r)
for r in reqRows:
    print("d=%.2g F=%.2f Fo=%.4g t=%.4g H=%.4g" % r)
