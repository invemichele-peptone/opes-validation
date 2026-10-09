"""Offline replay of the OPES kernel compression, to compare merge rules on the same deposited kernels.

The deposited kernels (KERNELS file) of OPES_METAD_EXPLORE do not depend on the bias, so replaying
them with different merge rules gives the bias each rule would produce along the same trajectory.
The compression mirrors OPESmetad::addKernel and OPESmetad::mergeKernels (recursive merge included).
"""

import numpy as np

# merge rules: (weights used for centre and variance, merged height)
#   PLAIN: heights, h1 + h2 (v2.9)
#   SPLIT: masses h*prod(sigma), mass / prod(merged sigma) (Sodkomkham et al.)
#   MASS:  heights for centre and variance, but the height conserves the mass
#   SHAPE: masses for centre and variance, but the height is h1 + h2
RULES = ("PLAIN", "SPLIT", "MASS", "SHAPE")


def wrap(d, periods):
    """Minimum-image difference, for the periodic CVs (period 0 means not periodic)."""
    out = np.array(d, dtype=float)
    for i, p in enumerate(periods):
        if p:
            out[..., i] -= p * np.round(out[..., i] / p)
    return out


class Compressor:
    def __init__(self, ncv, periods, threshold=1.0, rule="PLAIN", recursive=True, capacity=1000):
        assert rule in RULES
        self.ncv, self.periods, self.rule, self.recursive = ncv, np.asarray(periods, float), rule, recursive
        self.threshold2 = threshold**2
        self.c = np.empty((capacity, ncv))
        self.s = np.empty((capacity, ncv))
        self.h = np.empty(capacity)
        self.n = 0

    def _grow(self):
        for name in ("c", "s", "h"):
            a = getattr(self, name)
            setattr(self, name, np.concatenate([a, np.empty_like(a)]))

    def mergeable(self, center, exclude=-1):
        """Index of the closest kernel within the threshold, in units of its sigma, or -1."""
        if self.n == 0:
            return -1
        d = wrap(center - self.c[: self.n], self.periods) / self.s[: self.n]
        norm2 = np.sum(d * d, axis=1)
        if exclude >= 0:
            norm2[exclude] = np.inf
        k = int(np.argmin(norm2))  # first index of the minimum, as in PLUMED
        return k if norm2[k] < self.threshold2 else -1

    def merge(self, k1, c2, s2, h2):
        """Merge kernel (c2, s2, h2) into kernel k1."""
        c1, s1, h1 = self.c[k1].copy(), self.s[k1].copy(), self.h[k1]
        c1 = c2 + wrap(c1 - c2, self.periods)  # fix PBC, as in mergeKernels
        if self.rule in ("SPLIT", "SHAPE"):
            m1, m2 = h1 * np.prod(s1), h2 * np.prod(s2)
        else:
            m1, m2 = h1, h2
        m = m1 + m2
        w1, w2 = m1 / m, m2 / m
        delta = c2 - c1
        c = c1 + w2 * delta
        s = np.sqrt(w1 * s1**2 + w2 * s2**2 + w1 * w2 * delta**2)
        for i, p in enumerate(self.periods):
            if p:  # bring back in the periodic domain [-p/2, p/2)
                c[i] -= p * np.floor((c[i] + p / 2) / p)
        if self.rule == "PLAIN":
            h = h1 + h2
        elif self.rule == "SPLIT":
            h = m / np.prod(s)
        elif self.rule == "MASS":
            h = (h1 * np.prod(s1) + h2 * np.prod(s2)) / np.prod(s)
        else:  # SHAPE
            h = h1 + h2
        self.c[k1], self.s[k1], self.h[k1] = c, s, h

    def add(self, center, sigma, height):
        taker = self.mergeable(center)
        if taker < 0:
            if self.n == len(self.h):
                self._grow()
            self.c[self.n], self.s[self.n], self.h[self.n] = center, sigma, height
            self.n += 1
            return
        self.merge(taker, center, sigma, height)
        if not self.recursive:
            return
        giver = taker
        taker = self.mergeable(self.c[giver], exclude=giver)
        while taker >= 0:
            if taker > giver:
                taker, giver = giver, taker
            self.merge(taker, self.c[giver].copy(), self.s[giver].copy(), self.h[giver])
            for a in (self.c, self.s, self.h):  # erase giver
                a[giver : self.n - 1] = a[giver + 1 : self.n]
            self.n -= 1
            giver = taker
            taker = self.mergeable(self.c[giver], exclude=giver)

    @property
    def kernels(self):
        return self.c[: self.n], self.s[: self.n], self.h[: self.n]


def replay(deposited, ncv, periods, rule, threshold=1.0, until=None, snapshots=None):
    """Compress the deposited kernels (rows of a KERNELS file) with a merge rule.

    Returns the compressed kernels at time `until` (or at the end), or a list of them at the
    times in `snapshots`.
    """
    comp = Compressor(ncv, periods, threshold, rule)
    out, snaps = [], list(snapshots or [])
    for row in deposited:
        if until is not None and row[0] > until + 1e-9:
            break
        while snaps and row[0] > snaps[0] + 1e-9:
            out.append((snaps.pop(0), *[a.copy() for a in comp.kernels]))
        comp.add(row[1 : 1 + ncv], row[1 + ncv : 1 + 2 * ncv], row[1 + 2 * ncv])
    if snapshots is not None:
        out += [(t, *[a.copy() for a in comp.kernels]) for t in snaps]
        return out
    return comp.kernels


def kde(points, kernels, periods, cutoff, chunk=200_000):
    """Sum of truncated (energy-shifted) kernels at the points, as in OPES with CUTOFF_SCHEME=PLAIN or PR."""
    c, s, h = kernels
    val_at_cutoff = np.exp(-0.5 * cutoff**2)
    points = np.atleast_2d(points)
    out = np.empty(len(points))
    step = max(1, chunk // max(1, len(h)))
    for a in range(0, len(points), step):
        d = wrap(points[a : a + step, None, :] - c[None, :, :], periods) / s[None, :, :]
        norm2 = np.sum(d * d, axis=2)
        out[a : a + step] = np.sum(h * np.where(norm2 < cutoff**2, np.exp(-0.5 * norm2) - val_at_cutoff, 0), axis=1)
    return out


def bias(points, kernels, periods, cutoff, kdenorm, epsilon, prefactor, no_zed=False):
    """OPES bias a*log(P/Z+epsilon) at the points (a = prefactor, in energy units), and Z."""
    p = kde(points, kernels, periods, cutoff) / kdenorm
    zed = 1.0 if no_zed else np.mean(kde(kernels[0], kernels, periods, cutoff)) / kdenorm
    return prefactor * np.log(p / zed + epsilon), zed
