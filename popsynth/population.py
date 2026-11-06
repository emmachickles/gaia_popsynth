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
CO_KSTARS = (13, 14)  # neutron star, black hole


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


def evolve(ib, nproc=2, seed=1, **bse_overrides):
    """Evolve with COSMIC default flags plus overrides (e.g. alpha1=[3.0, 3.0]).

    `seed` fixes COSMIC's per-binary random draws (natal kicks, etc.), so a rerun
    with the same seed reproduces the population exactly.
    """
    bse = default_bse_dict(**bse_overrides)
    bpp, bcm, initC, kick = Evolve.evolve(initialbinarytable=ib, BSEDict=bse, nproc=nproc,
                                          randomseed=seed)
    return bpp, bcm, initC


def classify_final(bpp, bcm, initC):
    """Final state of every simulated binary, one row per bin_num.

    Categories, applied in this order:
      merged, disrupted        bin_state 1 / 2 (from COSMIC)
      interacting              still bound, a star fills its Roche lobe today
      wd_ms                    detached WD + main-sequence star
      wd_wd, wd_co             detached WD with another WD or a NS/BH companion
      wd_other                 detached WD with a giant, helium star, etc.
      ms_ms                    detached pair of main-sequence stars (not yet evolved)
      no_wd_other              detached, no WD: NS/BH pairs, giants, He stars
    """
    final = bcm.groupby("bin_num").tail(1).set_index("bin_num")
    k1, k2 = final.kstar_1, final.kstar_2
    is_wd1, is_wd2 = k1.isin(WD_KSTARS), k2.isin(WD_KSTARS)
    is_ms1, is_ms2 = k1.isin(MS_KSTARS), k2.isin(MS_KSTARS)
    is_co1, is_co2 = k1.isin(CO_KSTARS), k2.isin(CO_KSTARS)
    bound = final.bin_state == 0
    interacting = bound & ((final.RRLO_1 >= 1) | (final.RRLO_2 >= 1))
    detached = bound & ~interacting

    state = np.full(len(final), "no_wd_other", dtype=object)
    state[(detached & is_ms1 & is_ms2).values] = "ms_ms"
    state[(detached & ((is_wd1 & is_co2) | (is_wd2 & is_co1))).values] = "wd_co"
    state[(detached & is_wd1 & is_wd2).values] = "wd_wd"
    wd_ms = detached & ((is_wd1 & is_ms2) | (is_wd2 & is_ms1))
    state[wd_ms.values] = "wd_ms"
    wd_any = detached & (is_wd1 | is_wd2)
    state[(wd_any & ~wd_ms & ~(is_wd1 & is_wd2) & ~(is_co1 | is_co2)).values] = "wd_other"
    state[interacting.values] = "interacting"
    state[(final.bin_state == 2).values] = "disrupted"
    state[(final.bin_state == 1).values] = "merged"
    return pd.Series(state, index=final.index, name="final_state")


def final_table(bpp, bcm, initC):
    """All binaries, one row each: final state plus present-day properties."""
    final = bcm.groupby("bin_num").tail(1).set_index("bin_num")
    ic = initC.set_index("bin_num").loc[final.index]
    out = pd.DataFrame({
        "final_state": classify_final(bpp, bcm, initC),
        "kstar_1": final.kstar_1, "kstar_2": final.kstar_2,
        "mass_1": final.mass_1, "mass_2": final.mass_2,
        "porb": final.porb, "ecc": final.ecc, "tphys": final.tphys,
        "bin_state": final.bin_state, "merger_type": final.merger_type,
        "m1_zams": ic.mass_1, "m2_zams": ic.mass_2, "porb_zams": ic.porb,
    }, index=final.index)
    out["had_ce"] = out.index.isin(set(bpp.loc[bpp.evol_type == CE_EVOL_TYPE, "bin_num"]))
    return out


def select_wdms(bpp, bcm, initC):
    """Keep present-day detached WD + MS binaries; orient columns as MS / WD."""
    final = bcm.groupby("bin_num").tail(1).set_index("bin_num")
    state = classify_final(bpp, bcm, initC)
    keep = state == "wd_ms"
    f = final[keep]
    wd1 = f.kstar_1.isin(WD_KSTARS).values  # True where star 1 is the WD

    def pick(col):
        a, b = f[f"{col}_1"].values, f[f"{col}_2"].values
        return np.where(wd1, b, a), np.where(wd1, a, b)  # (ms, wd)

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


def build_all(n_binaries=100_000, seed=1, nproc=2, **bse_overrides):
    """Stage 1 with every binary kept. Returns (final_table, wdms_table, mass_formed)."""
    ib, mass_formed = sample_initial(n_binaries, seed=seed, nproc=nproc)
    bpp, bcm, initC = evolve(ib, nproc=nproc, seed=seed, **bse_overrides)
    final = final_table(bpp, bcm, initC)
    wdms = select_wdms(bpp, bcm, initC)
    wdms.attrs["mass_formed"] = mass_formed
    wdms.attrs["n_binaries"] = n_binaries
    wdms.attrs["bse_overrides"] = {k: str(v) for k, v in bse_overrides.items()}
    return final, wdms


def build(n_binaries=100_000, seed=1, nproc=2, **bse_overrides):
    """Stage 1: returns the WD+MS DataFrame only (see build_all for everything)."""
    return build_all(n_binaries, seed=seed, nproc=nproc, **bse_overrides)[1]
