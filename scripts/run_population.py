"""Stage 1 driver: build and cache a COSMIC WD+MS population.

Usage: python scripts/run_population.py --tag fiducial --alpha 1.0 --n 60000
"""
import argparse
import time

from popsynth import analysis, population

ap = argparse.ArgumentParser()
ap.add_argument("--tag", required=True)
ap.add_argument("--n", type=int, default=60000, help="number of zero-age binaries")
ap.add_argument("--alpha", type=float, default=1.0, help="common envelope efficiency alpha1")
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--nproc", type=int, default=2)
a = ap.parse_args()

t = time.time()
pop = population.build(a.n, seed=a.seed, nproc=a.nproc, alpha1=[a.alpha, a.alpha])
runtime = round(time.time() - t, 1)
analysis.save_population(pop, a.tag, alpha=a.alpha, seed=a.seed, runtime_s=runtime)
print(f"{a.tag}: {len(pop)} WD+MS systems from {a.n} binaries in {runtime} s -> data/pop_{a.tag}.csv.gz")
