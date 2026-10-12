"""Replays the OpenFOAM parcel update for a free 1D particle in still salt (Uc = 0):
    x  += U_old*dt                         (KinematicParcel::move: s = trackTime*U_ before calc)
    U  += (A - beta*U)*dtEff               (integrationScheme::delta, A = Su/m, beta = 1/tauP)
with the BoundedBrownianMotion amplitude (f = sqrt(m kT (1-exp(-2x)))/dtEff) or the stock
Li-Ahmadi one (f^2 = 2 m kT/(tauP dt)), and the analytical or CrankNicolson dtEff.
Prints the velocity std and the numerical diffusivity MSD/(2t) against D = kT tauP/m."""
import math
import numpy as np

kB, T, mu = 1.380649e-23, 946.15, 4.51e-3
rhop = 10980.0


def run(d, dt, nstep, amp, scheme, n=4000, seed=1):
    rng = np.random.default_rng(seed)
    m = rhop * math.pi / 6 * d ** 3
    tau = rhop * d ** 2 / (18 * mu)
    beta, x = 1 / tau, dt / tau
    dtEff = (-math.expm1(-beta * dt)) / beta if scheme == "analytical" else dt / (1 + beta * dt / 2)
    if amp == "bounded":
        om = -math.expm1(-2 * x) if x > 1e-12 else 2 * x
        f = math.sqrt(m * kB * T * om) / dtEff
    else:
        f = math.sqrt(2 * m * kB * T / (tau * dt))
    U = np.zeros(n)
    X = np.zeros(n)
    for _ in range(nstep):
        X += U * dt
        A = f * rng.standard_normal(n) / m
        U += (A - beta * U) * dtEff
    D = kB * T * tau / m
    t = nstep * dt
    return dict(x=x, Ustd=U.std(), vth=math.sqrt(kB * T / m), Dnum=(X ** 2).mean() / (2 * t), D=D)


if __name__ == "__main__":
    print(f"{'d [m]':>7} {'dt [s]':>7} {'dt/tauP':>9} {'scheme':>12} {'amplitude':>8} "
          f"{'U std':>9} {'v_th':>9} {'D_num/D':>10}")
    for d in (1e-9, 1e-7, 1e-5):
        for dt in (1e-6, 1e-3):
            for scheme in ("analytical", "CrankNicolson"):
                for amp in ("bounded", "stock"):
                    r = run(d, dt, 400, amp, scheme)
                    print(f"{d:7.0e} {dt:7.0e} {r['x']:9.2e} {scheme:>12} {amp:>8} "
                          f"{r['Ustd']:9.3g} {r['vth']:9.3g} {r['Dnum'] / r['D']:10.3g}")
