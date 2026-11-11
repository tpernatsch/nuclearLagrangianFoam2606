// V0: unit test of newmanSphere.H against an independent long-double series
// Build/run: g++ -O2 -I../../massTransferModels/gasDiffusionSimple Test-newmanSphere.C && ./a.out

#include "newmanSphere.H"
#include <cstdio>
#include <cmath>
#include <initializer_list>

using namespace newmanSphere;

// Series with many terms in long double (independent of the switch)
static long double remRef(long double Fo)
{
    const long double p = 3.14159265358979323846264338327950288L;
    long double s = 0;
    for (long n = 200000; n >= 1; --n)          // small terms first
    {
        s += expl(-(long double)n*n*p*p*Fo)/((long double)n*n);
    }
    return 6.0L/(p*p)*s;
}

static long double ShRef(long double Fo)
{
    const long double p = 3.14159265358979323846264338327950288L;
    long double s1 = 0, s2 = 0;
    for (long n = 200000; n >= 1; --n)
    {
        const long double e = expl(-(long double)n*n*p*p*Fo);
        s1 += e;
        s2 += e/((long double)n*n);
    }
    return (2.0L*p*p/3.0L)*s1/s2;
}

int main()
{
    int fails = 0;
    auto check = [&](const char* what, double Fo, double got, double ref, double tol)
    {
        const double err = std::fabs(got - ref)/std::fabs(ref);
        const bool ok = err < tol;
        if (!ok) ++fails;
        std::printf("%-10s Fo=%-9.3e got=%-.12e ref=%-.12e rel=%.2e %s\n",
            what, Fo, got, ref, err, ok ? "ok" : "FAIL");
    };

    // F and Sh across both branches, incl. both sides of the switch
    const double Fos[] =
        {1e-6, 1e-4, 1e-3, 0.01, 0.0499999, 0.05, 0.0500001, 0.1, 0.3, 1.0, 2.0};

    for (double Fo : Fos)
    {
        const double R = (double)remRef(Fo);
        check("F", Fo, F(Fo), 1.0 - R, 1e-6);
        check("1-F", Fo, remaining(Fo), R, 1e-6);
        check("Sh", Fo, Sh(Fo), (double)ShRef(Fo), 1e-6);
    }

    // Limits
    check("Sh0", 1e-10, Sh(1e-10), 2.0/std::sqrt(pi*1e-10), 1e-4);
    check("ShInf", 3.0, Sh(3.0), ShInf, 1e-9);

    // Step average reproduces the Newman loss: prod over steps of
    // exp(-1.5 ShMean dFo) = 1 - F, from Fo = 0 (singular first step)
    for (double FoEnd : {1e-4, 0.05, 0.3})
    {
        const int N = 1000;
        double lnR = 0;
        for (int i = 0; i < N; ++i)
        {
            const double a = FoEnd*i/N, b = FoEnd*(i + 1)/N;
            lnR -= 1.5*ShMean(a, b)*(b - a);
        }
        check("ShMean", FoEnd, std::exp(lnR), (double)remRef(FoEnd), 1e-9);
    }

    // Explicit-step Sh: prod over steps of (1 - 1.5 ShStep dFo) = 1 - F
    for (double FoEnd : {1e-4, 0.05, 0.3})
    {
        const int N = 300;
        double R = 1;
        for (int i = 0; i < N; ++i)
        {
            const double a = FoEnd*i/N, b = FoEnd*(i + 1)/N;
            R *= 1 - 1.5*ShStep(a, b)*(b - a);
        }
        check("ShStep", FoEnd, R, (double)remRef(FoEnd), 1e-9);
    }

    // First-step average is finite and equals the exact first-step loss
    const double Fo1 = 1e-8;
    check("ShMean0", Fo1, ShMean(0, Fo1),
        -std::log((double)remRef(Fo1))/(1.5*Fo1), 1e-6);

    std::printf("%s (%d failures)\n", fails ? "FAILED" : "PASSED", fails);
    return fails ? 1 : 0;
}
