#!/usr/bin/env python3
"""
claddingAnalysis.py - thermocouple profile T(z), cladding heat generation and the axial vs radial
conduction check in the Hastelloy cladding, for the heated-pool capsules.

    python3 claddingAnalysis.py              tables on screen + analysis/A_profile.png

Model of the cladding (steady, 1D radial, k constant per wall, uniform volumetric generation q_c):
    (1/r) d/dr (k r dT/dr) + q_c = 0,   T(r_o) = T_o,   -k dT/dr (r_i) = q_i  (flux leaving the liquid)
    T_w = T(r_i) = T_o + q_i r_i ln(r_o/r_i)/k + q_c/(4k) [r_o^2 - r_i^2 - 2 r_i^2 ln(r_o/r_i)]
i.e. q_i = h_c (T_w - T_o - dT_gen),  h_c = k/(r_i ln(r_o/r_i)),  dT_gen = q_c/(4k)[...]
Flat plate of thickness t (q_i into the plate at the liquid side, T_o on the far side):
    T_w = T_o + q_i t/k + q_c t^2/(2k)
"""
import math
import os

import numpy as np

from tprofile import profile

HERE = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------------------------------- data
TC_Z = [0.0, 37.5e-3, 75.0e-3]           # A: thermocouples on the axis inside the inner ring [m]
TC_T = [620.0, 690.0, 710.0]             # [C]

K_CLAD = 22.6                            # Hastelloy N [W/m/K] (~670 C)
RHO_CLAD = 8860.0                        # [kg/m3]
GEN_CLAD = 2.0                           # cladding heat generation [W/g]
Q_SALT = 1.0e8                           # salt heat generation [W/m3]

GEO = {
    "A": dict(H=75e-3, r_liq_in=2.5e-3, r_liq_out=4.3e-3, r_clad_out=8.5e-3, r_ring_in=1.25e-3, t_bottom=5e-3),
    "B": dict(H=35e-3, r_liq_in=0.0, r_liq_out=3.5e-3, r_clad_out=8.5e-3, r_ring_in=None, t_bottom=5e-3),
}


# ----------------------------------------------------------------------------- interpolants
def _vectorised(method):
    g = profile(TC_Z, TC_T, method)
    return np.vectorize(g), np.vectorize(g.d1), np.vectorize(g.d2)


PROFILES = {           # same profiles setupCase.py uses (tprofile.py)
    "natural cubic spline": _vectorised("natural"),
    "spline, dT/dz = 0 at the top": _vectorised("clamped_top"),
    "monotone (PCHIP)": _vectorised("pchip"),
}


# ----------------------------------------------------------------------------- cladding generation
def cladding_numbers(geo):
    g = GEO[geo]
    qc = GEN_CLAD * RHO_CLAD * 1e3                  # W/g -> W/m3
    H = g["H"]
    ri, ro = g["r_liq_out"], g["r_clad_out"]
    out = dict(qc=qc)
    out["V_salt"] = math.pi * (g["r_liq_out"] ** 2 - g["r_liq_in"] ** 2) * H
    out["Q_salt"] = Q_SALT * out["V_salt"]
    out["Q_outer"] = qc * math.pi * (ro ** 2 - ri ** 2) * H
    out["Q_bottom"] = qc * math.pi * ro ** 2 * g["t_bottom"]
    out["Q_ring"] = qc * math.pi * (g["r_liq_in"] ** 2 - g["r_ring_in"] ** 2) * H if g["r_ring_in"] else 0.0
    out["h_c"] = K_CLAD / (ri * math.log(ro / ri))
    out["dT_gen"] = qc / (4 * K_CLAD) * (ro ** 2 - ri ** 2 - 2 * ri ** 2 * math.log(ro / ri))
    out["dT_gen_plate"] = qc * g["t_bottom"] ** 2 / (2 * K_CLAD)
    # flux leaving the liquid through the outer wall if all salt + ring heat goes that way
    Aw = 2 * math.pi * ri * H
    out["q_i"] = (out["Q_salt"] + out["Q_ring"]) / Aw
    out["dT_clad"] = out["q_i"] / out["h_c"] + out["dT_gen"]
    out["q_o"] = (out["q_i"] * ri + qc * (ro ** 2 - ri ** 2) / 2) / ro    # at the outer surface
    if g["r_ring_in"]:
        out["q_ring_into_salt"] = out["Q_ring"] / (2 * math.pi * g["r_liq_in"] * H)
    return out


# ----------------------------------------------------------------------------- axial vs radial
def axial_ratio(d2f, z, area, radial_per_length):
    """Net axial conduction into a slice, k A |T''|, over the radial heat per unit length."""
    return K_CLAD * area * np.abs(d2f(z)) / radial_per_length


def main():
    a = cladding_numbers("A")
    b = cladding_numbers("B")
    print("Cladding heat generation (2 W/g x 8.86 g/cm3 = %.3g W/m3 = %.1f W/cm3)" % (a["qc"], a["qc"] / 1e6))
    print(f"{'':28}{'A':>10}{'B':>10}")
    for key, label, unit in [("Q_salt", "salt", "W"), ("Q_outer", "outer cladding (to 17 mm)", "W"),
                             ("Q_ring", "inner ring (A)", "W"), ("Q_bottom", "bottom plate (17 mm, 5 mm)", "W"),
                             ("h_c", "cladding conductance h_c", "W/m2K"),
                             ("q_i", "flux leaving the liquid", "W/m2"), ("q_o", "flux at cladding outer", "W/m2"),
                             ("dT_gen", "dT_gen (wall, own heat)", "K"), ("dT_clad", "T_w - T_o total", "K"),
                             ("dT_gen_plate", "dT_gen (bottom plate)", "K")]:
        print(f"  {label:<26}{a[key]:>10.4g}{b[key]:>10.4g}  {unit}")
    print(f"  inner ring heat INTO the salt (A): {a['q_ring_into_salt'] / 1e3:.3g} kW/m2 "
          f"(the ring is not adiabatic any more)")

    g = GEO["A"]
    z = np.linspace(0, g["H"], 301)
    area_outer = math.pi * (g["r_clad_out"] ** 2 - g["r_liq_out"] ** 2)
    area_ring = math.pi * (g["r_liq_in"] ** 2 - g["r_ring_in"] ** 2)
    radial_outer = (a["Q_salt"] + a["Q_ring"] + a["Q_outer"]) / g["H"]      # W/m leaving radially
    radial_ring = a["Q_ring"] / g["H"]
    print("\nA: axial / radial conduction (k A |d2T/dz2| over radial heat per length), using the profile")
    print(f"  {'profile':<30}{'outer clad max':>15}{'mean':>8}{'inner ring max':>16}{'mean':>8}")
    ratios = {}
    for name, (f, df, d2f) in PROFILES.items():
        ro = axial_ratio(d2f, z, area_outer, radial_outer)
        rr = axial_ratio(d2f, z, area_ring, radial_ring)
        ratios[name] = ro
        print(f"  {name:<30}{100 * ro.max():>14.2f}%{100 * ro.mean():>7.2f}%{100 * rr.max():>15.2f}%{100 * rr.mean():>7.2f}%")
    print("  (end values of a spline's 2nd derivative are set by its end condition, not by data)")

    plot(z, ratios)


def plot(z, ratios):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = ["#2a78d6", "#eb6834", "#1baf7a"]          # categorical slots 1-3 (light)
    ink, ink2, grid, surf = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": grid, "axes.labelcolor": ink2,
                         "xtick.color": ink2, "ytick.color": ink2, "text.color": ink})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 5.2), sharey=True, facecolor=surf,
                                   gridspec_kw=dict(width_ratios=[1.4, 1]))
    zz = 1e3 * z
    for (name, (f, _, _)), c in zip(PROFILES.items(), colors):
        ax1.plot(f(z), zz, color=c, lw=2, label=name)
    ax1.scatter(TC_T, [1e3 * v for v in TC_Z], s=46, color=ink, zorder=5, label="thermocouples",
                edgecolor=surf, linewidth=2)
    for t, zt in zip(TC_T, TC_Z):
        ax1.annotate(f"{t:.0f} °C", (t, 1e3 * zt), xytext=(8, -3 if zt < 0.07 else -10),
                     textcoords="offset points", color=ink2, fontsize=8)
    ax1.set_xlabel("temperature [°C]")
    ax1.set_ylabel("z [mm]  (A, inner ring / axis)")
    ax1.set_title("Thermocouple profile T(z)", loc="left", color=ink, fontsize=10)
    ax1.legend(frameon=False, loc="lower right", fontsize=8)

    for (name, r), c in zip(ratios.items(), colors):
        ax2.plot(100 * r, zz, color=c, lw=2)
    ax2.set_xlabel("axial / radial conduction [%]  (outer cladding)")
    ax2.set_title("Axial conduction check", loc="left", color=ink, fontsize=10)
    for ax in (ax1, ax2):
        ax.set_facecolor(surf)
        ax.grid(True, color=grid, lw=0.6)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    ax1.set_ylim(-2, 77)
    fig.tight_layout()
    os.makedirs(os.path.join(HERE, "analysis"), exist_ok=True)
    out = os.path.join(HERE, "analysis", "A_profile.png")
    fig.savefig(out, dpi=150, facecolor=surf)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
