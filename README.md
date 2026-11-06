# gaia_popsynth — prototype

Forward model from binary population synthesis to the Gaia astrometric orbit
catalog, starting with white dwarf + main-sequence (WD+MS) binaries. The goal
is a pipeline where each stage can be swapped out: a better Galaxy model,
gaiamock in place of the analytic selection, and DR4 in place of DR3.

```
COSMIC (BSE)            toy disk + dust           simplified Gaia             compare
Moe & Di Stefano  --->  positions, A_G, G  --->   a0, sigma(G), DR3 cuts  ---> P, M_WD, counts
constant SFR 10 Gyr     photocentre flux ratio    AMRF dark-companion triage   vs DR3 / published
popsynth/population.py  popsynth/galaxy.py        popsynth/gaia.py             scripts/compare.py
```

## Install

```bash
git clone https://github.com/<you>/gaia_popsynth && cd gaia_popsynth
pip install -e ".[notebooks,data]"
```

## Notebooks

The notebooks are the main walkthrough. They're committed with outputs, so they read on GitHub without running anything.

| Notebook | What it does |
| --- | --- |
| [01_cosmic_population](notebooks/01_cosmic_population.ipynb) | Stage 1: COSMIC WD+MS population, BSE flags, where orbits end up (ZAMS vs present-day period, CE vs no CE) |
| [02_galaxy_and_gaia_selection](notebooks/02_galaxy_and_gaia_selection.ipynb) | Stages 2–3: placing systems in the Galaxy, normalisation check, the DR3 selection funnel, AMRF triage |
| [03_compare_dr3_dr4](notebooks/03_compare_dr3_dr4.ipynb) | Stage 4: physics variants vs published DR3 results, and a DR4 prediction |

They use the two cached populations in `data/` (60,000 zero-age binaries each), so they run in about a minute.

## Command-line pipeline

```bash
python scripts/run_population.py --tag fiducial --alpha 1.0 --n 60000 --seed 11   # ~3 min on 2 cores
python scripts/run_population.py --tag alpha5   --alpha 5.0 --n 60000 --seed 11
python scripts/fetch_dr3_orbits.py      # optional: needs Gaia archive access
python scripts/compare.py fiducial alpha5   # figures/ + data/summary.json
```

Both variants use the same seed, so they start from identical zero-age binaries
and differ only in common envelope efficiency (alpha1).

## First results (60,000 zero-age binaries per variant)

| Quantity | alpha = 1 (COSMIC default) | alpha = 5 (efficient CE) | Published reference |
| --- | --- | --- | --- |
| Present-day WD+MS with P = 100–1000 d (fraction of all WD+MS) | 0.7% | 4.7% | — |
| Intrinsic dN/dP slope, 100–1000 d | −2.7 (71 sim systems) | −1.85 (502) | −0.25 (Yamaguchi+2025) |
| DR3 detections passing triage | ~5,700 | ~28,500 | ~3,200 (Shahaf+2024) |
| …of which P = 100–1000 d | 35% | 80% | nearly all |
| Median M_WD of DR3 detections | 0.59 | 0.64 | peak ~0.6 (Yamaguchi+2025) |
| DR4 detections passing triage | ~45,000 | ~87,000 | — |

What this already shows, with toy-model caveats:

1. **The default model has a gap where Gaia finds WD+MS orbits.** With alpha = 1, the
   100–1000 d range is nearly empty, which is the known problem with standard
   population synthesis. Raising alpha to 5 fills it.
2. **Efficient CE fills the window but gives the wrong shape.** The predicted period
   distribution falls far too steeply (−1.85 vs −0.25). This matches the conclusion
   of Yamaguchi+2025 that the alpha formalism does not reproduce the data, and
   points to AGB mass transfer stability rather than CE efficiency alone.
3. **DR4 is dominated by systems that never interacted.** DR4's longer baseline adds
   mostly P > 1000 d systems that never went through CE (CE fraction among
   detections drops from 31% to 7% at alpha = 1). That's a clean, testable
   DR4 prediction about the boundary between interacting and non-interacting systems.

**Do not trust absolute counts yet** (factor of several). The model produces
about 1.2 WD+MS per pc² versus roughly 2.7 white dwarfs of all kinds per pc² locally,
which is plausible to within about 2×. More important, the triage here (companion darker than any single MS
star) is looser than the Shahaf class-III selection plus color-excess cut used for the
observed ~3,200. Yamaguchi+2025 found only about 7% of mock WD+MS orbits survive the real triage.
Compare shapes, not totals, until the triage is matched. The toy period
efficiency also removes about 98% of systems that pass the signal cuts (most have P > 1000 d),
so DR3 counts depend heavily on it. Replacing it with gaiamock is the most important upgrade.

## What is simplified (and the planned replacement)

| Stage | Prototype | Replace with |
| --- | --- | --- |
| Initial binaries | COSMIC multidim (Moe & Di Stefano 2017), solar Z, constant SFR 10 Gyr | Metallicity spread, age–metallicity relation |
| Galaxy | Uniform cylinder R < 2 kpc, scale height 100→350 pc with age, Sigma_formed = 50 Msun/pc² | Galaxia or FIRE-2 m12i (ananke) + mwdust Combined19 |
| Photometry | Andrae+2018 BC_G(Teff) for both stars; WD Teff clipped to 8,000 K | WD cooling-model photometry (e.g. Bédard+2020), MIST for MS |
| Astrometric precision | sigma_fov = 0.1 mas × 10^(0.2(G−13)), N_fov = 40 (DR3) / 78 (DR4) | Holl+2023 relation; real scan times and angles (GOST) |
| Detection | Analytic a0/sigma, DR3 cuts, period efficiency 1/(1+(P/T)^6), 1-yr dip | **gaiamock** epoch simulation + DR3 cascade (El-Badry+2024) |
| Triage | AMRF above single-MS maximum, photometric mass from M_G | Shahaf+2019/2023 classes + color-excess cut |
| Statistics | Each simulated system cloned ~1,500×, so detected samples are lumpy | 10–100× more zero-age binaries (cluster), or importance sampling in P |

`popsynth.gaia.select(mock, release)` is the interface to keep. A gaiamock-based
version only has to return the same `detected` / `wd_class` columns.

## Next steps

1. Get observed DR3 orbits (`scripts/fetch_dr3_orbits.py`, run where the Gaia
   archive is reachable). They're pushed through the same triage code.
2. Swap in gaiamock for `gaia.select`; check the mock reproduces the DR3 period
   peak near 550 d before trusting any physics conclusion.
3. Vary the physics the literature argues about: AGB mass transfer stability
   (`qcflag` / `qcrit_array`), CE efficiency, WD natal kicks (eccentricities).
4. DR4: rerun with `GaiaRelease.dr4()` and the 66-month scanning law. Validate on
   the DR4 pre-release epoch astrometry sample.

## References

El-Badry et al. 2024 (gaiamock), arXiv:2411.00088 · Yamaguchi, El-Badry & Shahaf 2025,
arXiv:2505.14786 · Shahaf et al. 2024 · Halbwachs et al. 2023 · Shahaf et al. 2019 (AMRF) ·
Andrae et al. 2018 (BC_G) · Breivik et al. 2020 (COSMIC), arXiv:1911.00903 ·
Moe & Di Stefano 2017 · Hollands et al. 2018 (local WD density).
