"""Summarize build/time_nlist_cases.sh: cost of OPES with and without the exact neighbor list.

usage: python build/summarize_nlist_cases.py [--pool] <outdir> [<outdir> ...]
With several directories, each case is shown once per directory, or once in all with --pool.
"""
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

VARIANTS = ["none", "exact", "exact_p2", "exact_p5", "exact_p8"]


def load(dirs, pool):
    data = defaultdict(list)
    for d in dirs:
        for line in open(Path(d) / "times.dat"):
            if line.startswith("#"):
                continue
            c, variant, rep, nsteps, wall, opes, nker, nlker, nlsteps, *update = line.split()
            if len(dirs) > 1 and not pool:
                c = f"{c} ({Path(d).name})"
            data[(c, variant)].append(dict(nsteps=int(nsteps), wall=float(wall), opes=float(opes),
                                           nker=int(nker), nlker=float(nlker), nlsteps=float(nlsteps),
                                           update=float(update[0]) if update else np.nan))
    return data


def main():
    args = sys.argv[1:]
    pool = "--pool" in args
    data = load([a for a in args if a != "--pool"], pool)
    cases = sorted({c for c, _ in data}, key=lambda c: (c.split("-")[0][-2:], c))
    present = [v for v in VARIANTS if any(k[1] == v for k in data)]
    print(f"{'case':34s} {'kernels':>8s} {'listed':>7s} {'rebuild':>8s}   OPES us/step: "
          + " ".join(f"{v:>15s}" for v in present) + f"   {'none/exact':>10s}   {'wall: none':>10s} {'exact':>7s} {'ratio':>6s}"
          + f"   {'update: none':>12s} {'exact':>6s}")
    for c in cases:
        runs = {v: data.get((c, v), []) for v in VARIANTS}
        if not runs["none"] or not runs["exact"]:
            continue
        cost = {v: np.array([r["opes"] / r["nsteps"] * 1e6 for r in rs]) for v, rs in runs.items() if rs}
        wall = {v: np.array([r["wall"] for r in rs]) for v, rs in runs.items() if rs}
        update = {v: np.mean([r["update"] / r["nsteps"] * 1e6 for r in rs]) for v, rs in runs.items() if rs}
        nker = np.mean([r["nker"] for r in runs["none"]])
        listed = np.mean([r["nlker"] for r in runs["exact"]])
        rebuild = np.mean([r["nlsteps"] for r in runs["exact"]])
        cells = [f"{cost[v].mean():7.3f} ± {cost[v].std(ddof=1):5.3f}" if v in cost else f"{'-':>15s}" for v in present]
        print(f"{c:34s} {nker:8.0f} {listed:7.0f} {rebuild:8.0f}                 " + " ".join(cells) + "   "
              f"{cost['none'].mean() / cost['exact'].mean():10.2f}   "
              f"{wall['none'].mean():9.1f}s {wall['exact'].mean():6.1f}s {wall['none'].mean() / wall['exact'].mean():6.2f}"
              f"   {update['none']:12.3f} {update['exact']:6.3f}")


if __name__ == "__main__":
    main()
