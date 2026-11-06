"""Stage 3: simplified Gaia astrometric-orbit selection (DR3-like, DR4-ready).

This is an analytic stand-in for gaiamock (El-Badry et al. 2024), which
simulates epoch astrometry along the real scanning law and refits it with the
DR3 cascade. Swap `select` for a gaiamock-based function with the same
signature when that code is available; DR4 is the same call with
GaiaRelease.dr4().

What is modelled
----------------
* Photocentre semi-major axis a0 = a_rel * |m_WD/(m_WD+m_MS) - F_WD/(F_WD+F_MS)|.
* Per-transit along-scan precision sigma_fov(G) (rough fit to the DR3 trend).
* Orbit significance a0/sigma_a0 with sigma_a0 ~ sigma_fov * sqrt(2/N_fov).
* DR3 orbital-solution cuts (Halbwachs et al. 2023), in simplified form:
    a0/sigma_a0 > 158/sqrt(P/d),  a0/sigma_a0 > 5,
    parallax/sigma_parallax > 20000 d / P,  G < 19.
* Period coverage: efficiency falls off for P beyond ~ the baseline.
* WD+MS triage via the astrometric mass ratio function (AMRF, Shahaf+2019):
    the companion must be too dark to be a main-sequence star of that mass.

What is NOT modelled (and matters): scan-angle/time sampling, eccentricity
dependence, acceleration-solution competition, spurious solutions, RUWE
pre-selection, triples. Expect counts good to a factor of ~2, shapes better.
"""
from dataclasses import dataclass

import numpy as np

G_NEWTON_AU = 39.478  # AU^3 / (Msun yr^2), i.e. 4 pi^2
RSUN_AU = 0.00465047


@dataclass
class GaiaRelease:
    name: str
    baseline_d: float  # time baseline in days
    n_fov: int  # typical number of field-of-view transits
    g_lim: float = 19.0

    @classmethod
    def dr3(cls):
        return cls("DR3", baseline_d=1038.0, n_fov=40)  # 34 months

    @classmethod
    def dr4(cls):
        return cls("DR4", baseline_d=2006.0, n_fov=78)  # 66 months


def sigma_fov(g):
    """Rough per-FoV-transit along-scan precision [mas] versus G.

    ~0.1 mas floor for G < 13, roughly x1.6 per magnitude fainter
    (~1.6 mas at G = 19). Replace with the Holl+2023 relation used in gaiamock.
    """
    g = np.asarray(g, float)
    return 0.1 * 10 ** (0.2 * np.clip(g - 13.0, 0, None))


def photocentre_a0(mock):
    """Photocentre semi-major axis [mas] (angular)."""
    a_rel_au = mock.sep.values * RSUN_AU
    f = mock.m_wd / (mock.m_wd + mock.m_ms)
    b = mock.flux_ratio / (1 + mock.flux_ratio)
    return a_rel_au * np.abs(f - b) * mock.parallax.values  # mas


def period_efficiency(p, baseline):
    """Smooth fall-off of orbit recoverability with period (toy).

    1 / (1 + (P/baseline)^6): ~1 well inside the baseline, 0.5 at P = baseline,
    ~0.1 at 1.5x baseline (DR3 has almost no orbits beyond ~1500 d), plus a
    deficit near 1 yr, where orbit and parallax are degenerate.
    """
    p = np.asarray(p, float)
    eff = 1.0 / (1.0 + (p / baseline) ** 6)
    yr = 365.25
    eff *= 1 - 0.6 * np.exp(-0.5 * ((p - yr) / 15.0) ** 2)
    return eff


def amrf(a0, parallax, m_ms, p_d):
    """Astrometric mass-ratio function A = (a0/plx) / (M1^(1/3) P_yr^(2/3))."""
    return (a0 / parallax) / (m_ms ** (1 / 3) * (p_d / 365.25) ** (2 / 3))


def mass_from_abs_g(m_g):
    """Crude main-sequence mass from absolute G: 1 Msun at M_G = 4.7, 0.6 at ~7.7.

    Used identically on the mock and on real DR3 data, so the triage is a
    like-for-like forward model even though the relation itself is rough.
    """
    return np.clip(10 ** (-0.075 * (np.asarray(m_g, float) - 4.7)), 0.1, 3.0)


def amrf_ms_max(m_ms, lum_exp=4.0):
    """Largest AMRF any single MS companion can give (G-band flux ~ M^lum_exp)."""
    q = np.linspace(0.01, 1.0, 200)[None, :]
    s = q**lum_exp
    a = q / (1 + q) ** (2 / 3) * np.abs(1 - s * (1 + q) / (q * (1 + s)))
    return np.full(np.shape(m_ms), a.max())


def select(mock, release=None, rng=None):
    """Return mock with detection columns and a boolean 'detected'."""
    release = release or GaiaRelease.dr3()
    rng = np.random.default_rng(rng)
    m = mock.copy()
    a0 = photocentre_a0(m)
    sig = sigma_fov(m.G.values)
    sig_a0 = sig * np.sqrt(2.0 / release.n_fov)
    sig_plx = 1.5 * sig / np.sqrt(release.n_fov)
    snr = a0 / sig_a0
    p = m.porb.values

    m["a0"] = a0
    m["a0_snr"] = snr
    m["plx_snr"] = m.parallax.values / sig_plx
    m["m1_phot"] = mass_from_abs_g(m.M_G.values)
    m["amrf"] = amrf(a0, m.parallax.values, m.m1_phot.values, p)

    cuts = (
        (m.G.values < release.g_lim)
        & (snr > 5)
        & (snr > 158.0 / np.sqrt(p))
        & (m.plx_snr.values > 20000.0 / p)
    )
    eff = period_efficiency(p, release.baseline_d)
    detected = cuts & (rng.uniform(size=len(m)) < eff)
    m["detected"] = detected
    # triage: dark companion (AMRF above the single-MS-companion maximum)
    m["wd_class"] = detected & (m.amrf.values > amrf_ms_max(m.m1_phot.values))
    m.attrs.update(mock.attrs)
    m.attrs["release"] = release.name
    return m
