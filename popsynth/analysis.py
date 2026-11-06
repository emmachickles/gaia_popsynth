"""Shared helpers for loading cached populations and summarising mock catalogs."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from . import galaxy, gaia

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA = REPO_ROOT / "data"

# Published reference points (not fitted).
REF = {
    "n_wdms_dr3": 3200,  # ~3,200 WD + luminous-companion astrometric binaries (Shahaf+2024)
    "dndp_slope": -0.25,  # intrinsic dN/dP ~ P^-0.25 on 100-1000 d (Yamaguchi+2025)
    "m_wd_peak": 0.6,  # Msun (Yamaguchi+2025)
}

COLORS = {"fiducial": "#2a6fdb", "alpha5": "#d1495b"}
LABELS = {"fiducial": "COSMIC default (alpha = 1)", "alpha5": "efficient CE (alpha = 5)"}


def save_population(pop, tag, **meta):
    DATA.mkdir(exist_ok=True)
    pop.to_csv(DATA / f"pop_{tag}.csv.gz")
    (DATA / f"pop_{tag}.json").write_text(json.dumps(dict(pop.attrs, tag=tag, **meta), indent=2))


def load_population(tag):
    pop = pd.read_csv(DATA / f"pop_{tag}.csv.gz", index_col=0)
    pop.attrs.update(json.loads((DATA / f"pop_{tag}.json").read_text()))
    return pop


def load_observed():
    """Observed DR3 dark-companion orbits, if fetch_dr3_orbits.py has been run."""
    path = DATA / "dr3_astrometric_orbits.csv.gz"
    if not path.exists():
        return None
    obs = pd.read_csv(path)
    return obs[obs.wd_class]


def run_pipeline(tag, seed=7):
    """Stages 2-3 for a cached population: realise in the Galaxy, select for DR3 and DR4."""
    pop = load_population(tag)
    mock = galaxy.realise(pop, rng=seed)
    out = {"pop": pop, "mock": mock}
    for rel in (gaia.GaiaRelease.dr3(), gaia.GaiaRelease.dr4()):
        out[rel.name] = gaia.select(mock, rel, rng=seed + 1)
    return out


def slope_100_1000(p, lo=100.0, hi=1000.0):
    """Max-likelihood power-law index of dN/dP between lo and hi days."""
    p = np.asarray(p)
    p = p[(p > lo) & (p < hi)]
    if len(p) < 5:
        return np.nan, len(p)
    grid = np.linspace(-3, 2, 501)
    ll = []
    for g in grid:
        norm = np.log(hi / lo) if abs(g + 1) < 1e-9 else (hi ** (g + 1) - lo ** (g + 1)) / (g + 1)
        ll.append(np.sum(g * np.log(p)) - len(p) * np.log(norm))
    return float(grid[int(np.argmax(ll))]), int(len(p))


def summarise(result):
    """Headline numbers for one variant's pipeline output."""
    pop, w = result["pop"], result["mock"].attrs["weight"]
    in_win = (pop.porb > 100) & (pop.porb < 1000)
    s = {
        "weight": w,
        "n_sim_wdms": len(pop),
        "mass_formed_sim": pop.attrs["mass_formed"],
        "frac_wdms_sim_P_100_1000": float(in_win.mean()),
        "n_sim_P_100_1000": int(in_win.sum()),
        "frac_had_ce": float(pop.had_ce.mean()),
    }
    for rel in ("DR3", "DR4"):
        d = result[rel][result[rel].wd_class]
        s[f"n_detected_{rel}"] = int(len(d))
        s[f"n_unique_sim_{rel}"] = int(d.sim_id.nunique())
        if len(d):
            s[f"P_median_{rel}"] = float(np.median(d.porb))
            s[f"frac_P_100_1000_{rel}"] = float(((d.porb > 100) & (d.porb < 1000)).mean())
            s[f"frac_had_ce_{rel}"] = float(d.had_ce.mean())
            s[f"m_wd_median_{rel}"] = float(np.median(d.m_wd))
            s[f"dist_median_{rel}"] = float(np.median(d.dist))
    s["intrinsic_dndp_slope_100_1000"], s["intrinsic_n_100_1000"] = slope_100_1000(pop.porb.values)
    return s
