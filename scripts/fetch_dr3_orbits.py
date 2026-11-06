"""Fetch Gaia DR3 astrometric orbits and apply the SAME triage as the mock.

Run this on a machine that can reach the Gaia archive (the Claude sandbox
cannot). Requires astroquery.

    python scripts/fetch_dr3_orbits.py          # writes data/dr3_astrometric_orbits.csv.gz

Then scripts/compare.py picks the file up automatically.
"""
import numpy as np
from astroquery.gaia import Gaia

from popsynth.analysis import DATA
from popsynth.gaia import amrf, amrf_ms_max, mass_from_abs_g

QUERY = """
SELECT n.source_id, n.nss_solution_type, n.period, n.period_error, n.eccentricity,
       n.a_thiele_innes, n.b_thiele_innes, n.f_thiele_innes, n.g_thiele_innes,
       n.parallax AS plx_orbit, n.parallax_error AS plx_orbit_error, n.goodness_of_fit,
       s.phot_g_mean_mag, s.bp_rp, s.ruwe, s.l, s.b
FROM gaiadr3.nss_two_body_orbit AS n
JOIN gaiadr3.gaia_source AS s USING (source_id)
WHERE n.nss_solution_type IN ('Orbital', 'AstroSpectroSB1')
"""


def thiele_innes_a0(A, B, F, G):
    """Photocentre semi-major axis from Thiele-Innes constants (Halbwachs+2023)."""
    u = (A**2 + B**2 + F**2 + G**2) / 2
    v = A * G - B * F
    return np.sqrt(u + np.sqrt((u + v) * (u - v)))


def main():
    Gaia.ROW_LIMIT = -1
    t = Gaia.launch_job_async(QUERY).get_results().to_pandas()
    t["a0"] = thiele_innes_a0(t.a_thiele_innes, t.b_thiele_innes, t.f_thiele_innes, t.g_thiele_innes)
    t["M_G"] = t.phot_g_mean_mag + 5 * np.log10(t.plx_orbit / 100.0)  # no extinction correction
    t["m1_phot"] = mass_from_abs_g(t.M_G)
    t["amrf"] = amrf(t.a0, t.plx_orbit, t.m1_phot, t.period)
    t["wd_class"] = t.amrf > amrf_ms_max(t.m1_phot)
    out = DATA / "dr3_astrometric_orbits.csv.gz"
    t.to_csv(out, index=False)
    print(f"{len(t)} orbits, {t.wd_class.sum()} pass the dark-companion triage -> {out}")


if __name__ == "__main__":
    main()
