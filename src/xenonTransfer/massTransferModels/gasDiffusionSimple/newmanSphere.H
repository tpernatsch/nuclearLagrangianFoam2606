/*---------------------------------------------------------------------------*\
    Newman (1931) diffusion out of a sphere, zero surface concentration,
    uniform initial content (Crank, The Mathematics of Diffusion, eqs 6.20-6.21)

    Fo = 4*D*t/d^2 (radius-based Fo = D*t/R^2)
    F(Fo)  = 1 - (6/pi^2) sum exp(-n^2 pi^2 Fo)/n^2       removed fraction
           = 6 sqrt(Fo/pi) - 3 Fo  + O(exp(-1/Fo))        short time
    Bulk-referred Sh = k_L d/D with dYb/dt = -(6/d) k_L Yb:
        dln(1 - F)/dFo = -1.5 Sh
    Sh -> 2/sqrt(pi Fo) (Fo -> 0), 2 pi^2/3 (Fo -> inf)

    Header-only, plain double, so it can be unit-tested outside OpenFOAM.
\*---------------------------------------------------------------------------*/

#ifndef newmanSphere_H
#define newmanSphere_H

#include <cmath>

namespace newmanSphere
{

constexpr double pi = 3.14159265358979323846;

// Short-time branch below this Fo (relative error in F ~1e-10 at the switch)
constexpr double FoSwitch = 0.05;

// Long-time Sherwood number 2 pi^2/3
constexpr double ShInf = 2.0*pi*pi/3.0;

// Remaining fraction 1 - F
inline double remaining(const double Fo)
{
    if (Fo <= 0) return 1.0;

    if (Fo < FoSwitch)
    {
        return 1.0 - 6.0*std::sqrt(Fo/pi) + 3.0*Fo;
    }

    double s = 0;
    for (int n = 1; n <= 30; ++n)
    {
        s += std::exp(-n*n*pi*pi*Fo)/(n*n);
    }
    return 6.0/(pi*pi)*s;
}

// Removed fraction F
inline double F(const double Fo)
{
    if (Fo <= 0) return 0.0;

    if (Fo < FoSwitch)
    {
        return 6.0*std::sqrt(Fo/pi) - 3.0*Fo;
    }

    return 1.0 - remaining(Fo);
}

// Instantaneous bulk-referred Sherwood number
inline double Sh(const double Fo)
{
    if (Fo <= 0) return HUGE_VAL;

    if (Fo < FoSwitch)
    {
        return (2.0/std::sqrt(pi*Fo) - 2.0)/remaining(Fo);
    }

    double s1 = 0, s2 = 0;
    for (int n = 1; n <= 30; ++n)
    {
        const double e = std::exp(-n*n*pi*pi*Fo);
        s1 += e;
        s2 += e/(n*n);
    }
    return ShInf*s1/s2;
}

// Sherwood number averaged over [Fo0, Fo1]: the exact Newman loss of the
// interval, ln((1-F0)/(1-F1))/(1.5 dFo); finite at Fo0 = 0
inline double ShMean(const double Fo0, const double Fo1)
{
    const double a = Fo0 > 0 ? Fo0 : 0;
    const double b = Fo1;

    if (b <= a) return Sh(b);

    // Cancellation guard for tiny intervals
    if (b - a < 1e-8*b) return Sh(0.5*(a + b));

    return (std::log(remaining(a)) - std::log(remaining(b)))/(1.5*(b - a));
}

// Sherwood number for an explicit (forward Euler) source over [Fo0, Fo1]:
// 1 - 1.5 Sh dFo = (1-F1)/(1-F0), so the explicit step loses exactly the
// Newman fraction of the step
inline double ShStep(const double Fo0, const double Fo1)
{
    const double a = Fo0 > 0 ? Fo0 : 0;
    const double b = Fo1;

    if (b <= a) return Sh(b);
    if (b - a < 1e-8*b) return Sh(0.5*(a + b));

    return -std::expm1(std::log(remaining(b)) - std::log(remaining(a)))
        /(1.5*(b - a));
}

} // namespace newmanSphere

#endif
