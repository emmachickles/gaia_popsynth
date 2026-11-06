"""Stages 2-4 from the command line: place cached populations, select, summarise, plot.

    python scripts/compare.py fiducial alpha5

Writes figures/*.png and data/summary.json. Notebook 03 does the same step by step.
If data/dr3_astrometric_orbits.csv.gz exists (from fetch_dr3_orbits.py), the observed
dark-companion sample is overplotted.
"""
import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from popsynth import analysis, plots


def main(tags):
    res = {t: analysis.run_pipeline(t) for t in tags}
    obs = analysis.load_observed()
    summary = {t: analysis.summarise(r) for t, r in res.items()}
    summary["reference"] = analysis.REF
    (analysis.DATA / "summary.json").write_text(json.dumps(summary, indent=2, default=float))

    figs = analysis.REPO_ROOT / "figures"
    figs.mkdir(exist_ok=True)
    plots.period_distributions(res, obs).savefig(figs / "period_distributions.png", dpi=160)
    plots.period_vs_wd_mass(res).savefig(figs / "period_vs_wd_mass.png", dpi=160)
    plt.close("all")
    print(json.dumps(summary, indent=2, default=float))


if __name__ == "__main__":
    main(sys.argv[1:] or ["fiducial", "alpha5"])
