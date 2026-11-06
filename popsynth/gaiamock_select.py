"""Stage 3 with gaiamock: epoch astrometry along the real scanning law + DR3 cascade.

Drop-in alternative to `popsynth.gaia.select`. For each mock system it uses
gaiamock (El-Badry et al. 2024, https://github.com/kareemelbadry/gaiamock) to
generate per-transit along-scan astrometry at the real Gaia scan times and
angles, then fits it with the DR3 cascade (5-parameter -> acceleration ->
12-parameter orbit) and applies the DR3 orbital-solution cuts.

Setup (once):
    git clone https://github.com/kareemelbadry/gaiamock      # e.g. next to this repo
    cd gaiamock && gcc -shared -o kepler_solve_astrometry.so kepler_solve_astrometry.c \
        -lgsl -lgslcblas -lm -fPIC
    # healpix_scans/ (from the gaiamock README) must sit inside gaiamock/
    export GAIAMOCK_PATH=/path/to/gaiamock

Cost control
------------
gaiamock fits each source individually (~0.1-1 s), so it is run on a random
subsample of systems that pass a loose analytic pre-cut, and each fitted
system carries a weight (N passing pre-cut / N fitted). The pre-cut keeps
anything with analytic a0/sigma_a0 > 3, below the cascade's own threshold
of 5, so it removes almost nothing gaiamock would have accepted.

Sky sampling
------------
gaiamock looks up scan times per HEALPix (nside 64) cell. Positions can be
snapped to the centres of a coarser grid (`snap_nside`, default 8, i.e. 768
cells, ~7 deg), so only that subset of the 49,152 scan files is needed.
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from . import gaia

Z_SUN_GAIAMOCK = 20.8  # pc; gaiamock's convention for the Sun's height


def _import_gaiamock(path=None):
    path = Path(path or os.environ.get("GAIAMOCK_PATH", Path(__file__).resolve().parents[1] / "gaiamock"))
    if not (path / "gaiamock.py").exists():
        raise ImportError(f"gaiamock not found at {path}; set GAIAMOCK_PATH (see module docstring)")
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))
    import gaiamock  # noqa: E402

    return gaiamock, path


def scan_cells_needed(snap_nside=8):
    """nside-64 HEALPix indices of the scan files needed for a given snap grid."""
    import healpy as hp

    theta, phi = hp.pix2ang(snap_nside, np.arange(hp.nside2npix(snap_nside)))
    return np.unique(hp.ang2pix(64, theta, phi))


def snap_radec(ra, dec, snap_nside=8):
    """Move each position to the centre of its nside=snap_nside HEALPix cell."""
    import healpy as hp

    pix = hp.ang2pix(snap_nside, np.radians(90.0 - dec), np.radians(ra))
    theta, phi = hp.pix2ang(snap_nside, pix)
    return np.degrees(phi), 90.0 - np.degrees(theta)


def precut(mock, release="dr3", p_max=None, snr_min=3.0):
    """Loose analytic pre-selection before running the expensive cascade."""
    gm, _ = _import_gaiamock()
    n_fov = 40 if release == "dr3" else 78
    p_max = p_max or (3000.0 if release == "dr3" else 6000.0)
    sig_t = gm.al_uncertainty_per_ccd_interp(mock.G.values) / np.sqrt(8)
    snr = gaia.photocentre_a0(mock) / (sig_t * np.sqrt(2.0 / n_fov))
    return (mock.G.values < 19) & (mock.porb.values > 10) & (mock.porb.values < p_max) & (snr > snr_min)


def _worker(args):
    """Run gaiamock's cascade for one system (in a joblib worker)."""
    path, kw = args
    gm, _ = _import_gaiamock(path)
    np.random.seed(kw.pop("seed"))
    try:
        return gm.run_full_astrometric_cascade(c_funcs=None, verbose=False, show_residuals=False, **kw)
    except Exception:  # a failed fit counts as not detected
        return 23 * [0]


COLS = ["fit_plx", "fit_sig_plx", "A", "sig_A", "B", "sig_B", "F", "sig_F", "G_ti", "sig_G_ti",
        "fit_period", "fit_sig_period", "phi_p", "sig_phi_p", "fit_ecc", "fit_sig_ecc", "fit_inc",
        "fit_a0", "fit_sig_a0", "n_vis", "n_obs", "F2", "ruwe"]


def select(mock, release=None, rng=None, n_max=3000, n_jobs=2, snap_nside=8, gaiamock_path=None):
    """gaiamock version of gaia.select: returns the fitted subsample with weights.

    Output has the input columns plus the cascade results, `solution`
    ('orbit', 'accel7', 'accel9', '5par', 'none'), `detected` (DR3 orbital
    cuts), `wd_class` (AMRF dark-companion triage on the *fitted* orbit), and
    `w` (real systems represented by each row).
    """
    from joblib import Parallel, delayed

    rel = (release.name if hasattr(release, "name") else (release or "DR3")).lower()
    gm, path = _import_gaiamock(gaiamock_path)
    rng = np.random.default_rng(rng)

    keep = precut(mock, rel)
    cand = mock[keep]
    n_pre = len(cand)
    if n_pre == 0:
        raise ValueError("no systems pass the pre-cut")
    sub = cand.iloc[rng.choice(n_pre, size=min(n_max, n_pre), replace=False)].copy()
    sub["w"] = n_pre / len(sub)

    ra, dec = gm.xyz_to_radec(sub.x.values, sub.y.values, sub.z.values + Z_SUN_GAIAMOCK)
    if snap_nside:
        ra, dec = snap_radec(ra, dec, snap_nside)
    sub["ra"], sub["dec"] = ra, dec

    n = len(sub)
    jobs = []
    for i, r in enumerate(sub.itertuples()):
        jobs.append((str(path), dict(
            ra=float(ra[i]), dec=float(dec[i]), parallax=float(r.parallax), pmra=0.0, pmdec=0.0,
            m1=float(r.m_ms), m2=float(r.m_wd), period=float(r.porb),
            Tp=float(rng.uniform(0, r.porb)), ecc=float(r.ecc),
            omega=float(rng.uniform(0, 2 * np.pi)), inc_deg=float(np.degrees(np.arccos(r.cos_i))),
            w=float(rng.uniform(0, 2 * np.pi)), phot_g_mean_mag=float(r.G), f=float(r.flux_ratio),
            data_release=rel, seed=int(rng.integers(2**31)),
        )))
    res = Parallel(n_jobs=n_jobs)(delayed(_worker)(j) for j in jobs)
    fit = pd.DataFrame(np.array(res, dtype=float), columns=COLS, index=sub.index)
    sub = pd.concat([sub, fit], axis=1)

    flag = sub.fit_plx.values
    sub["solution"] = np.select([flag == 0, flag == -1, flag == -9, flag == -7, flag > 0],
                                ["none", "5par", "accel9", "accel7", "orbit"], "orbit")
    orb = sub.solution == "orbit"
    P = np.where(orb, sub.fit_period, np.nan)
    a0_snr = sub.fit_a0 / sub.fit_sig_a0
    plx_snr = sub.fit_plx / sub.fit_sig_plx
    # DR3 orbital-solution publication cuts (Halbwachs+2023; as in gaiamock)
    det = orb & (a0_snr > 5) & (a0_snr > 158 / np.sqrt(P)) & (plx_snr > 20000 / P) \
        & (sub.fit_sig_ecc < 0.079 * np.log(P) - 0.244) & (sub.F2 < 25)
    sub["detected"] = det.fillna(False).astype(bool)

    m1_phot = gaia.mass_from_abs_g(sub.M_G.values)
    amrf = np.where(orb, gaia.amrf(sub.fit_a0.values, np.where(orb, sub.fit_plx, 1.0), m1_phot,
                                   np.where(orb, sub.fit_period, 1.0)), np.nan)
    sub["m1_phot"], sub["amrf"] = m1_phot, amrf
    sub["wd_class"] = sub.detected & (sub.amrf > gaia.amrf_ms_max(m1_phot))
    sub.attrs.update(mock.attrs)
    sub.attrs.update(release=rel.upper(), selector="gaiamock", n_precut=n_pre, n_fitted=n)
    return sub
