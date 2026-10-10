"""Compare the outputs of build/replay_nlist.sh.

usage: python build/check_replay.py <replay-dir> [<other-replay-dir>]

With one directory: the largest difference between each neighbor-list variant (nlp: v2.9,
nla: EXACT rebuild, nlab: EXACT rebuild + FULL zed) and the reference without a list.
With two: the largest difference between the same action in the two directories, e.g. to
check an MPI or OpenMP replay against a serial one.
Bias and force differences are in kBT (T=300 K, kJ/mol); zed differences are relative.
"""
import sys
from pathlib import Path

import numpy as np

KBT = 0.0083144626 * 300
ACTIONS = ["ref", "nlp", "nla", "nlab"]


def load(d, name):
    bias = np.loadtxt(d / f"bias{name}.out")
    der = np.loadtxt(d / f"der{name}.out")
    nder = int(der[:, 1].max()) + 1
    forces = der[:, 2:6].reshape(-1, nder, 4)  # frame, cv, action
    return bias, forces


def diffs(b1, f1, b2, f2, i, j):
    """largest |Δbias|, |Δforce| (kBT) and relative |Δzed| of action i in run 1 vs action j in run 2"""
    n = min(len(b1), len(b2))
    db = np.abs(b1[:n, 1 + i] - b2[:n, 1 + j]).max() / KBT
    df = np.abs(f1[:n, :, i] - f2[:n, :, j]).max() / KBT
    dz = np.abs(b1[:n, 5 + i] / b2[:n, 5 + j] - 1).max()
    return db, df, dz


def outputs(d):
    names = sorted({p.name[len("bias"):-len(".out")] for p in d.glob("bias*.out")})
    return [(d, n) for n in names]


def main():
    dirs = [Path(a) for a in sys.argv[1:]]
    print(f"{'case':28s} " + "  ".join(f"{a:>28s}" for a in (ACTIONS[1:] if len(dirs) == 1 else ACTIONS)))
    print(f"{'':28s} " + "  ".join(f"{'bias   force    zed':>28s}" for _ in (ACTIONS[1:] if len(dirs) == 1 else ACTIONS)))
    for case in sorted(p for p in dirs[0].iterdir() if p.is_dir()):
        for d, name in outputs(case):
            b1, f1 = load(d, name)
            if len(dirs) == 1:
                cells = [diffs(b1, f1, b1, f1, j, 0) for j in range(1, 4)]
            else:
                b2, f2 = load(dirs[1] / case.name, name)
                cells = [diffs(b1, f1, b2, f2, j, j) for j in range(4)]
            print(f"{case.name + name:28s} " + "  ".join(f"{db:8.1e} {df:8.1e} {dz:8.1e}" for db, df, dz in cells))


if __name__ == "__main__":
    main()
