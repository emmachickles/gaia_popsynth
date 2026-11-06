"""Stage 1: COSMIC population of present-day WD + MS binaries.

Draws zero-age binaries from the Moe & Di Stefano (2017) multidimensional
sampler, gives them a constant star formation history, evolves them with
COSMIC/BSE, and keeps systems that are a detached white dwarf + main-sequence
pair today.

Output columns (one row per simulated system):
    m_ms, m_wd [Msun], porb [d], sep [Rsun], ecc,
    lum_ms, teff_ms, lum_wd, teff_wd, t_cool [Myr], age [Myr],
    m1_zams, m2_zams, porb_zams [d], had_ce (bool)
plus the population normalisation in DataFrame.attrs["mass_formed"] [Msun].
"""
import numpy as np
import pandas as pd
from cosmic.sample.initialbinarytable import InitialBinaryTable
from cosmic.evolve import Evolve

from .bse_defaults import default_bse_dict

WD_KSTARS = (10, 11, 12)
MS_KSTARS = (0, 1)
CE_EVOL_TYPE = 7  # bpp evol_type for common envelope


def sample_initial(n_binaries, seed=1, sf_duration=10_000.0, met=0.014, nproc=2):
    """Moe & Di Stefano initial binaries with constant SFR over sf_duration Myr.

    Returns (InitialBinaryTable, total stellar mass formed [Msun]).
    The mass formed includes the single stars COSMIC had to draw to produce
    this many binaries, so it is the right normalisation for absolute counts.
    """
    ib, m_singles, m_bin, n_s, n_b = InitialBinaryTable.sampler(
        "multidim",
        final_kstar1=list(WD_KSTARS),
        final_kstar2=list(MS_KSTARS),
        rand_seed=seed,
        nproc=nproc,
        SF_start=sf_duration,
        SF_duration=sf_duration,
        met=met,
        size=n_binaries,
        keep_singles=False,
    )
    return ib, float(m_singles + m_bin)


def evolve(ib, nproc=2, **bse_overrides):
    """Evolve with COSMIC default flags plus overrides (e.g. alpha1=[3.0, 3.0])."""
    bse = default_bse_dict(**bse_overrides)
    bpp, bcm, initC, kick = Evolve.evolve(initialbinarytable=ib, BSEDict=bse, nproc=nproc)
    return bpp, bcm, initC


def select_wdms(bpp, bcm, initC):
    """Keep present-day detached WD + MS binaries; orient columns as MS / WD."""
    final = bcm.groupby("bin_num").tail(1).set_index("bin_num")
    bound = final.bin_state == 0
    k1, k2 = final.kstar_1, final.kstar_2
    wd1 = k1.isin(WD_KSTARS) & k2.isin(MS_KSTARS)
    wd2 = k2.isin(WD_KSTARS) & k1.isin(MS_KSTARS)
    detached = (final.RRLO_1 < 1) & (final.RRLO_2 < 1)
    keep = bound & detached & (wd1 | wd2)
    f = final[keep]
    w1 = wd1[keep].values  # True where star 1 is the WD

    def pick(col):
        a, b = f[f"{col}_1"].values, f[f"{col}_2"].values
        return np.where(w1, b, a), np.where(w1, a, b)  # (ms, wd)

    m_ms, m_wd = pick("mass")
    lum_ms, lum_wd = pick("lum")
    teff_ms, teff_wd = pick("teff")
    ep_ms, ep_wd = pick("epoch")

    ic = initC.set_index("bin_num").loc[f.index]
    ce_bins = set(bpp.loc[bpp.evol_type == CE_EVOL_TYPE, "bin_num"])

    out = pd.DataFrame(
        {
            "m_ms": m_ms,
            "m_wd": m_wd,
            "porb": f.porb.values,
            "sep": f.sep.values,
            "ecc": f.ecc.values,
            "lum_ms": lum_ms,
            "teff_ms": teff_ms,
            "lum_wd": lum_wd,
            "teff_wd": teff_wd,
            # BSE epoch: time the current stellar type began, so t - epoch = WD cooling age
            "t_cool": f.tphys.values - ep_wd,
            "age": f.tphys.values,
            "m1_zams": ic.mass_1.values,
            "m2_zams": ic.mass_2.values,
            "porb_zams": ic.porb.values,
            "had_ce": np.isin(f.index.values, list(ce_bins)),
        },
        index=f.index,
    )
    return out


def build(n_binaries=100_000, seed=1, nproc=2, **bse_overrides):
    """Full stage 1: sample, evolve, select. Returns WD+MS DataFrame."""
    ib, mass_formed = sample_initial(n_binaries, seed=seed, nproc=nproc)
    bpp, bcm, initC = evolve(ib, nproc=nproc, **bse_overrides)
    pop = select_wdms(bpp, bcm, initC)
    pop.attrs["mass_formed"] = mass_formed
    pop.attrs["n_binaries"] = n_binaries
    pop.attrs["bse_overrides"] = {k: str(v) for k, v in bse_overrides.items()}
    return pop
