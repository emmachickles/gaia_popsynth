"""Stage 2: place simulated systems in a toy solar-neighbourhood disk.

Deliberately simple and swappable (target replacement: Galaxia, or FIRE-2
m12i via ananke, with mwdust Combined19 extinction).

Model
-----
* Star formation: constant over the population's SF duration, with surface
  density of mass ever formed SIGMA_FORMED at the Sun's radius.
* Geometry: uniform in a cylinder of radius R_MAX around the Sun; vertical
  sech^2-like (here exponential) profile whose scale height grows with age.
* Extinction: a dust layer with scale height H_DUST and midplane A_V per kpc.
* Photometry: Gaia G from bolometric luminosity and the Andrae+2018 BC_G(Teff)
  polynomial, for both components; the WD flux enters the photocentre.

Each simulated system is cloned n times (Poisson) so that the realisation is
an absolute-number mock catalog of the volume.
"""
import numpy as np
import pandas as pd

SIGMA_FORMED = 50.0  # Msun / pc^2 of stars ever formed near the Sun (assumption)
R_MAX = 2000.0  # pc, cylinder radius around the Sun
H_DUST = 100.0  # pc
AV_PER_KPC = 0.7  # mag / kpc in the midplane
AG_OVER_AV = 0.789  # approximate for G band
MBOL_SUN = 4.74


def scale_height(age_myr):
    """Age-dependent vertical scale height [pc]: ~100 pc young to ~350 pc at 10 Gyr."""
    return 100.0 + 25.0 * (np.asarray(age_myr) / 1000.0)


def bc_g(teff):
    """Gaia G bolometric correction, Andrae et al. 2018 (DR2 FLAME) polynomial.

    Valid 3300-8000 K; Teff is clipped to that range, which is crude for hot
    WDs but WDs barely move the photocentre in this regime.
    """
    t = np.clip(np.asarray(teff, dtype=float), 3300.0, 8000.0)
    d = t - 5772.0
    hot = np.array([6.000e-02, 6.731e-05, -6.647e-08, 2.859e-11, -7.197e-15])
    cool = np.array([1.749e00, 1.977e-03, 3.737e-07, -8.966e-11, -4.183e-14])
    c = np.where(t[..., None] >= 4000.0, hot, cool)
    return sum(c[..., i] * d**i for i in range(5))


def abs_g(lum, teff):
    lum = np.clip(np.asarray(lum, dtype=float), 1e-8, None)
    return MBOL_SUN - 2.5 * np.log10(lum) - bc_g(teff)


def extinction_g(d_pc, z_pc):
    """A_G integrated along a straight line from the Sun (z=0) to height z."""
    d = np.asarray(d_pc, float)
    z = np.abs(np.asarray(z_pc, float))
    sinb = np.where(d > 0, z / np.maximum(d, 1e-9), 0.0)
    # integral of exp(-s sinb / H) ds from 0 to d
    small = sinb * d / H_DUST < 1e-6
    path = np.where(small, d, H_DUST / np.maximum(sinb, 1e-12) * (1 - np.exp(-sinb * d / H_DUST)))
    return AG_OVER_AV * AV_PER_KPC * path / 1000.0


def realise(pop, rng=None, r_max=R_MAX, sigma_formed=SIGMA_FORMED, p_range=(10.0, 1e4)):
    """Clone systems into an absolute-number mock of the cylinder volume.

    Only systems with p_range[0] < P [d] < p_range[1] are cloned: outside that
    range Gaia cannot fit an astrometric orbit, and cloning them would only
    cost memory. Counts inside the window are unaffected.
    """
    rng = np.random.default_rng(rng)
    mass_real = sigma_formed * np.pi * r_max**2
    weight = mass_real / pop.attrs["mass_formed"]  # real systems per simulated one
    attrs = dict(pop.attrs)
    pop = pop[(pop.porb > p_range[0]) & (pop.porb < p_range[1])]
    pop.attrs.update(attrs)
    n = rng.poisson(weight, size=len(pop))
    mock = pop.loc[pop.index.repeat(n)].rename_axis("sim_id").reset_index()
    k = len(mock)

    r = r_max * np.sqrt(rng.uniform(size=k))
    phi = rng.uniform(0, 2 * np.pi, size=k)
    hz = scale_height(mock.age.values)
    z = rng.choice([-1, 1], size=k) * rng.exponential(hz)
    x, y = r * np.cos(phi), r * np.sin(phi)
    d = np.sqrt(x**2 + y**2 + z**2)

    mock["x"], mock["y"], mock["z"], mock["dist"] = x, y, z, d
    mock["parallax"] = 1000.0 / d  # mas

    ag = extinction_g(d, z)
    dm = 5 * np.log10(d) - 5
    g_ms = abs_g(mock.lum_ms, mock.teff_ms) + dm + ag
    g_wd = abs_g(mock.lum_wd, mock.teff_wd) + dm + ag
    f_ms, f_wd = 10 ** (-0.4 * g_ms), 10 ** (-0.4 * g_wd)
    mock["G"] = -2.5 * np.log10(f_ms + f_wd)
    mock["G_ms"] = g_ms
    mock["flux_ratio"] = f_wd / f_ms  # WD / MS in G
    # "observed" absolute G, NOT extinction-corrected, to match how
    # fetch_dr3_orbits.py computes it from real data
    mock["M_G"] = mock.G - dm
    mock["A_G"] = ag

    # orientation: isotropic inclination
    mock["cos_i"] = rng.uniform(-1, 1, size=k)
    mock.attrs.update(pop.attrs)
    mock.attrs["weight"] = weight
    mock.attrs["volume"] = dict(r_max=r_max, sigma_formed=sigma_formed)
    return mock
