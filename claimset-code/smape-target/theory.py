"""What the derivations predict, computed from each distribution exactly (by quadrature on its density, never from a
sample).

For positive outcomes, everything is computed on a fine grid of x = ln y. Expected sMAPE of the forecast f = e^c is
S(c) = 200 - 400 * integral of g(x) / (1 + e^{|x - c|}) dx  (since tanh(|z|/2) = 1 - 2 / (1 + e^{|z|})),
a convolution of the density g of ln y with a fixed kernel, computed for every c of the grid at once by FFT.
"""
import math

import numpy as np

N_GRID = 2 ** 19
K_FLAT = 800 * math.log(2)          # integral of 200 (1 - tanh(|z|/2)) dz over the real line (H14)


def _trapz_weights(n, dx):
    w = np.full(n, dx)
    w[0] = w[-1] = dx / 2
    return w


class Grid:
    """The density of ln y of a positive distribution on N_GRID points of its support."""

    def __init__(self, d, n=N_GRID):
        lo, hi = d.support
        self.x = np.linspace(lo, hi, n)
        self.dx = self.x[1] - self.x[0]
        g = d.pdf_x(self.x)
        self.w = _trapz_weights(n, self.dx)
        self.mass = float(g @ self.w)                 # 1 up to quadrature and truncation error
        self.g = g / self.mass
        self.gw = self.g * self.w
        self.cdf = np.concatenate([[0.0], np.cumsum((self.g[1:] + self.g[:-1]) / 2 * self.dx)])
        self.cdf /= self.cdf[-1]

    def quantile_x(self, levels):
        return np.interp(levels, self.cdf, self.x)

    def level_x(self, c):
        return np.interp(c, self.x, self.cdf)

    def density_x(self, c):
        return np.interp(c, self.x, self.g)

    def smape_all(self):
        """Expected sMAPE of e^c for every c of the grid, by FFT convolution."""
        n = self.x.size
        lags = np.arange(-(n - 1), n) * self.dx
        e = np.exp(-np.abs(lags))
        kern = e / (1.0 + e)                          # 1 / (1 + e^{|z|}), without overflow
        m = 1 << int(math.ceil(math.log2(3 * n)))
        conv = np.fft.irfft(np.fft.rfft(self.gw, m) * np.fft.rfft(kern, m), m)
        # conv[j + n - 1] = sum_i gw_i kern(x_j - x_i)
        return 200.0 - 400.0 * conv[n - 1: 2 * n - 1]

    def smape_at(self, c):
        """Expected sMAPE of e^c for given c, directly."""
        c = np.atleast_1d(np.asarray(c, dtype=float))
        out = np.empty(c.size)
        for i in range(0, c.size, 8):
            out[i:i + 8] = 200.0 * (np.tanh(np.abs(self.x[None, :] - c[i:i + 8, None]) / 2) @ self.gw)
        return out


def _golden(fn, a, b, tol=1e-10):
    phi = (math.sqrt(5) - 1) / 2
    c, d = b - phi * (b - a), a + phi * (b - a)
    fc, fd = fn(c), fn(d)
    while b - a > tol:
        if fc <= fd:
            b, d, fd = d, c, fc
            c = b - phi * (b - a)
            fc = fn(c)
        else:
            a, c, fc = c, d, fd
            d = a + phi * (b - a)
            fd = fn(d)
    return (a + b) / 2


def local_minima(xs, s, refine, rel=1e-9):
    """Interior local minima of s over the grid xs, each refined by golden-section search with `refine` between its
    neighbours; returns [(x, value)] in increasing order of value. Flat stretches (equal values within rel) count once."""
    out = []
    j = 1
    n = s.size
    while j < n - 1:
        if s[j] < s[j - 1]:
            k = j
            while k + 1 < n and abs(s[k + 1] - s[j]) <= rel * max(1.0, abs(s[j])):
                k += 1
            if k + 1 < n and s[k + 1] > s[j]:
                a, b = xs[j - 1], xs[k + 1]
                xm = _golden(lambda t: float(refine(t)), a, b, tol=(b - a) * 1e-9)
                out.append((xm, float(refine(xm))))
            j = k + 1
        else:
            j += 1
    return sorted(out, key=lambda t: t[1])


GL_X, GL_W = np.polynomial.legendre.leggauss(16)
PANELS = 3000


def _gauss(fn, a, b, panels):
    """Composite 16-point Gauss-Legendre integral of fn over [a, b]."""
    if b <= a:
        return 0.0
    edges = np.linspace(a, b, panels + 1)
    mid, half = (edges[1:] + edges[:-1]) / 2, (edges[1:] - edges[:-1]) / 2
    x = (mid[:, None] + half[:, None] * GL_X[None, :]).ravel()
    return float(fn(x) @ (half[:, None] * GL_W[None, :]).ravel())


class PositiveTheory:
    """Derived quantities of a positive or zero-inflated distribution."""

    def __init__(self, d):
        self.d = d
        self.p0 = d.p0
        self.grid = gr = Grid(d)
        self.mass = gr.mass                           # the density's integral over the support: 1 up to truncation
        self.S_plus = gr.smape_all()                  # expected sMAPE on the positive part, at every grid point
        # only the central part of the grid is searched: there the FFT result is accurate and minima can lie
        inner = (gr.cdf > 1e-6) & (gr.cdf < 1 - 1e-6)
        xs, sp = gr.x[inner], self.S_plus[inner]
        step = max(1, xs.size // 20000)
        self.minima_plus = local_minima(xs[::step], sp[::step], self.smape_exact)
        self.m_plus = self.minima_plus[0][1]
        self.median_x = float(gr.quantile_x(0.5))
        # moments: closed forms; the quadrature values are kept to check them (test_smape.py)
        self.mean_plus = d.moments.get("mean")
        self.second_plus = d.moments.get("second")
        self.inv_plus = d.moments.get("inv")

    def smape_exact(self, c):
        """Expected sMAPE on the positive part at the forecast e^c, by Gauss-Legendre quadrature of the density itself,
        split at c and at the density's kinks (more accurate than the grid where the minimum is flat)."""
        d, mass = self.d, self.grid.mass
        cuts = sorted({*d.support, *(k for k in d.kinks if d.support[0] < k < d.support[1]),
                       min(max(c, d.support[0]), d.support[1])})
        fn = lambda x: np.tanh(np.abs(x - c) / 2) * d.pdf_x(x)
        return 200.0 * sum(_gauss(fn, a, b, PANELS) for a, b in zip(cuts[:-1], cuts[1:])) / mass

    # ---------------------------------------------------------------- the target of sMAPE (H4-H8, H13)
    def smape_minima(self):
        """[(forecast, expected sMAPE)] of every local minimum over forecasts >= 0, lowest first. With zeros, the
        forecast 0 scores 200 (1 - p0) (a zero outcome forecast exactly scores 0, every other outcome 200)."""
        out = [(math.exp(c), self.p0 * 200 + (1 - self.p0) * v) for c, v in self.minima_plus]
        if self.p0 > 0:
            out.append((0.0, 200 * (1 - self.p0)))
        return sorted(out, key=lambda t: t[1])

    def zero_threshold(self):
        """H8: the share of zeros above which the forecast 0 has the lowest expected sMAPE."""
        return (200 - self.m_plus) / (400 - self.m_plus)

    def one_step(self):
        """H6: D (pull from above minus pull from below at the median) and the one-step quantile level
        1/2 + D / B'(1/2), with B'(1/2) = 1 - E[w(Z) tanh(|Z|/2)] / g(m) (z = x - ln m, w = sech^2(z/2) / 2)."""
        gr, c = self.grid, self.median_x
        z = gr.x - c
        w = 0.5 / np.cosh(z / 2) ** 2
        above, below = float((w * (z > 0)) @ gr.gw), float((w * (z < 0)) @ gr.gw)
        D = above - below
        slope = 1 - float((w * np.tanh(np.abs(z) / 2)) @ gr.gw) / float(gr.density_x(c))
        return D, 0.5 + D / slope

    def flat_approximation(self, c):
        """H14: 200 - 800 ln 2 * (density of ln y at c), on the positive part."""
        return 200 - K_FLAT * self.grid.density_x(c)

    def asymptotic_variances(self):
        """H11: large-sample n * Var of ln(sMAPE minimiser) at the global minimum, and n * Var of ln(sample median)."""
        gr = self.grid
        c = self.minima_plus[0][0]
        z = gr.x - c
        e_psi2 = 1e4 * float((1 / np.cosh(z / 2) ** 4) @ gr.gw)
        s2 = 200 * float(gr.density_x(c)) - 100 * float((np.tanh(np.abs(z) / 2) / np.cosh(z / 2) ** 2) @ gr.gw)
        v_smape = e_psi2 / s2 ** 2 if s2 > 0 else math.inf
        v_median = 1 / (4 * float(gr.density_x(self.median_x)) ** 2)
        return v_smape, v_median, s2

    # ---------------------------------------------------------------- expected scores of fixed forecasts
    def _lower(self, f, fn):
        """integral over positive y < f of fn(y) dF(y), for each f."""
        gr = self.grid
        out = np.empty(f.size)
        for i, ff in enumerate(f):
            if ff <= 0:
                out[i] = 0.0
                continue
            m = gr.x < math.log(ff)
            out[i] = float(fn(np.exp(gr.x[m]), ff) @ gr.gw[m])
        return out

    def expected(self, f):
        """Expected scores of forecasts f >= 0 under each measure (inf where infinite, H12), and finite ranking keys
        (differences in expected score are finite even where the scores are not): {measure: (scores, keys or None)}."""
        f = np.asarray(f, dtype=float)
        q = 1 - self.p0
        out = {}
        # RMSE
        if self.mean_plus is not None:
            mean = q * self.mean_plus
            key = (mean - f) ** 2
            if self.second_plus is not None:
                var = q * self.second_plus - mean ** 2
                out["RMSE"] = (np.sqrt(var + key), key)
            else:
                out["RMSE"] = (np.full(f.size, np.inf), key)
        else:
            out["RMSE"] = (np.full(f.size, np.inf), None)           # E Y infinite: no forecast is preferred
        # MAE: E|Y - f| = E Y - f + 2 E[(f - Y) 1{Y < f}]; key without E Y
        low = q * self._lower(f, lambda y, ff: ff - y) + self.p0 * f
        key = -f + 2 * low
        out["MAE"] = (key + (q * self.mean_plus if self.mean_plus is not None else np.inf), key)
        # MAPE
        if self.p0 > 0:
            out["MAPE"] = (np.full(f.size, np.inf), None)           # a zero outcome: infinite or undefined everywhere
        elif self.inv_plus is not None:
            key = 100 * (-f * self.inv_plus + 2 * self._lower(f, lambda y, ff: ff / y - 1))
            out["MAPE"] = (100 + key, key)
        else:                                                       # E[1/Y] infinite: only f = 0 is finite (100)
            sc = np.where(f == 0, 100.0, np.inf)
            out["MAPE"] = (sc, sc)
        # sMAPE
        sp = np.array([self.smape_exact(math.log(ff)) if ff > 0 else 200.0 for ff in f])
        sc = np.where(f > 0, 200 * self.p0 + q * sp, 200 * q)
        out["sMAPE"] = (sc, sc)
        return out

    def quantiles(self, levels):
        """Quantiles of the whole distribution (zeros included)."""
        levels = np.asarray(levels, dtype=float)
        if self.p0 == 0:
            return np.exp(self.grid.quantile_x(levels))
        inner = np.clip((levels - self.p0) / (1 - self.p0), 0, 1)
        return np.where(levels <= self.p0, 0.0, np.exp(self.grid.quantile_x(inner)))

    def level(self, f):
        """The share of the distribution at or below f."""
        f = np.asarray(f, dtype=float)
        pos = np.where(f > 0, self.grid.level_x(np.log(np.where(f > 0, f, 1.0))), 0.0)
        return self.p0 + (1 - self.p0) * pos

    def mape_minimiser(self):
        """The forecast minimising expected MAPE (H12), or None when MAPE targets nothing; 0 when E[1/Y] is infinite."""
        if self.p0 > 0:
            return None
        if self.inv_plus is None:
            return 0.0
        # minimiser: the median of the 1/y-weighted distribution
        gr = self.grid
        wi = gr.gw * np.exp(-gr.x)
        c = np.cumsum(wi) / wi.sum()
        return float(math.exp(np.interp(0.5, c, gr.x)))


class StraddlingTheory:
    """Normal outcomes of both signs: everything by quadrature in y."""

    def __init__(self, d, n=400001):
        self.d = d
        self.p0 = 0.0
        lo, hi = d.support
        self.y = np.linspace(lo, hi, n)
        self.w = _trapz_weights(n, self.y[1] - self.y[0])
        self.gw = d.pdf_y(self.y) * self.w
        self.mass = float(self.gw.sum())
        self.gw /= self.mass
        self.cdf = np.cumsum(self.gw)

    def smape_at(self, f):
        f = np.atleast_1d(np.asarray(f, dtype=float))
        out = np.empty(f.size)
        for i in range(0, f.size, 8):
            ff = f[i:i + 8, None]
            num, den = np.abs(self.y[None, :] - ff), np.abs(self.y[None, :]) + np.abs(ff)
            out[i:i + 8] = 200 * (np.where(num == 0, 0.0, num / np.where(den == 0, 1.0, den)) @ self.gw)
        return out

    def smape_minima(self):
        lo, hi = self.quantiles([0.001, 0.999])
        fs = np.unique(np.concatenate([np.linspace(lo, hi, 4001), [0.0]]))
        s = self.smape_at(fs)
        mins = local_minima(fs, s, lambda t: self.smape_at(t)[0])
        # the forecast 0 scores exactly 200 on every non-zero outcome; it is a candidate, never a local minimum here
        return mins

    def quantiles(self, levels):
        return np.interp(np.asarray(levels, dtype=float), self.cdf, self.y)

    def level(self, f):
        return np.interp(np.asarray(f, dtype=float), self.y, self.cdf)

    def expected(self, f):
        f = np.asarray(f, dtype=float)
        mu, sd = self.d.mu, self.d.sd
        key = (mu - f) ** 2
        mae = np.array([float(np.abs(self.y - ff) @ self.gw) for ff in f])
        mape = np.where(f == 0, 100.0, np.inf)             # E[1/|Y|] is infinite; the forecast 0 scores 100 exactly
        sm = self.smape_at(f)
        return {"RMSE": (np.sqrt(sd * sd + key), key), "MAE": (mae, mae), "MAPE": (mape, mape), "sMAPE": (sm, sm)}

    def mape_minimiser(self):
        return 0.0


def theory(d):
    return StraddlingTheory(d) if d.kind == "straddling" else PositiveTheory(d)
