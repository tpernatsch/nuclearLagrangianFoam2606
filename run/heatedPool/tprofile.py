"""
tprofile.py - temperature profiles T(z) through a few thermocouple readings (standard library only).

    f = profile(z, T, method)      f(z), f.d1(z), f.d2(z)

methods:
    "natural"     cubic spline, T'' = 0 at both ends
    "clamped_top" cubic spline, T'' = 0 at the bottom, dT/dz = 0 at the top (adiabatic free surface)
    "pchip"       monotone piecewise cubic (Fritsch-Carlson): never overshoots the readings
Outside the readings the end value is held (no extrapolation).
"""
import bisect


class _Cubic:
    """Piecewise cubic Hermite form: values y and end-of-interval slopes m on knots x."""

    def __init__(self, x, y, m):
        self.x, self.y, self.m = list(x), list(y), list(m)

    def _seg(self, z):
        x = self.x
        z = min(max(z, x[0]), x[-1])
        i = min(max(bisect.bisect_right(x, z) - 1, 0), len(x) - 2)
        h = x[i + 1] - x[i]
        return i, (z - x[i]) / h, h

    def __call__(self, z):
        i, t, h = self._seg(z)
        y, m = self.y, self.m
        return ((2 * t ** 3 - 3 * t ** 2 + 1) * y[i] + (t ** 3 - 2 * t ** 2 + t) * h * m[i]
                + (-2 * t ** 3 + 3 * t ** 2) * y[i + 1] + (t ** 3 - t ** 2) * h * m[i + 1])

    def d1(self, z):
        i, t, h = self._seg(z)
        y, m = self.y, self.m
        return ((6 * t ** 2 - 6 * t) * y[i] + (3 * t ** 2 - 4 * t + 1) * h * m[i]
                + (-6 * t ** 2 + 6 * t) * y[i + 1] + (3 * t ** 2 - 2 * t) * h * m[i + 1]) / h

    def d2(self, z):
        i, t, h = self._seg(z)
        y, m = self.y, self.m
        return ((12 * t - 6) * y[i] + (6 * t - 4) * h * m[i] + (-12 * t + 6) * y[i + 1]
                + (6 * t - 2) * h * m[i + 1]) / h ** 2


def _spline_slopes(x, y, top_slope=None):
    """Knot slopes of a cubic spline: natural at the bottom, natural or clamped (top_slope) at the top.
    Tridiagonal system in the slopes (C2 continuity), solved with the Thomas algorithm."""
    n = len(x)
    h = [x[i + 1] - x[i] for i in range(n - 1)]
    d = [(y[i + 1] - y[i]) / h[i] for i in range(n - 1)]
    a, b, c, r = [0.0] * n, [0.0] * n, [0.0] * n, [0.0] * n
    b[0], c[0], r[0] = 2.0, 1.0, 3.0 * d[0]                          # natural: T''(x0) = 0
    for i in range(1, n - 1):
        a[i], b[i], c[i] = h[i], 2.0 * (h[i - 1] + h[i]), h[i - 1]
        r[i] = 3.0 * (h[i] * d[i - 1] + h[i - 1] * d[i])
    if top_slope is None:                                           # natural: T''(xn) = 0
        a[-1], b[-1], r[-1] = 1.0, 2.0, 3.0 * d[-1]
    else:                                                           # clamped
        b[-1], r[-1] = 1.0, top_slope
    for i in range(1, n):                                           # Thomas
        w = a[i] / b[i - 1]
        b[i] -= w * c[i - 1]
        r[i] -= w * r[i - 1]
    m = [0.0] * n
    m[-1] = r[-1] / b[-1]
    for i in range(n - 2, -1, -1):
        m[i] = (r[i] - c[i] * m[i + 1]) / b[i]
    return m


def _pchip_slopes(x, y):
    n = len(x)
    h = [x[i + 1] - x[i] for i in range(n - 1)]
    d = [(y[i + 1] - y[i]) / h[i] for i in range(n - 1)]
    m = [0.0] * n
    for i in range(1, n - 1):
        if d[i - 1] * d[i] > 0:
            w1, w2 = 2 * h[i] + h[i - 1], h[i] + 2 * h[i - 1]
            m[i] = (w1 + w2) / (w1 / d[i - 1] + w2 / d[i])
    if n == 2:
        return [d[0], d[0]]
    for end in (0, -1):
        h0, h1 = (h[0], h[1]) if end == 0 else (h[-1], h[-2])
        d0, d1 = (d[0], d[1]) if end == 0 else (d[-1], d[-2])
        s = ((2 * h0 + h1) * d0 - h0 * d1) / (h0 + h1)
        if s * d0 <= 0:
            s = 0.0
        elif d0 * d1 < 0 and abs(s) > abs(3 * d0):
            s = 3 * d0
        m[end] = s
    return m


def profile(z, T, method="clamped_top"):
    pts = sorted(zip(z, T))
    x, y = [p[0] for p in pts], [p[1] for p in pts]
    if len(x) < 2:
        raise ValueError("a profile needs at least 2 readings")
    if method == "natural":
        m = _spline_slopes(x, y)
    elif method == "clamped_top":
        m = _spline_slopes(x, y, top_slope=0.0)
    elif method == "pchip":
        m = _pchip_slopes(x, y)
    else:
        raise ValueError(f"unknown profile method {method}")
    return _Cubic(x, y, m)
