"""MD tests of OPES with MPI and multiple walkers (WALKERS_MPI), with GROMACS built with MPI.

usage: python build/walkers_mpi.py run   OUTDIR LABEL [LABEL ...]   # MPI builds of PLUMED in opt/plumed
       python build/walkers_mpi.py check OUTDIR LABEL [LABEL ...]

Alanine dipeptide in vacuum, OPES_METAD and OPES_METAD_EXPLORE inputs of the Explore paper,
with GROMACS (gmx_mpi, opt/gromacs/<version>-mpi) running several walkers with -multidir:
- layout "walkers": 4 walkers sharing the bias (WALKERS_MPI), 1 MPI rank each, 5 ns per walker,
  SETS independent sets of walkers (replicas 4s..4s+3 of systems/ala2/tpr);
- layout "dd": 2 walkers sharing the bias, each with 2 MPI ranks (domain decomposition)
  and PLUMED_NUM_THREADS=2, 1 ns per walker, one set.
Each layout runs with the neighbor list (as in the inputs) and without it ("nolist": NLIST removed,
and NLIST_OFF for PLUMED >= 2.11, where the list is used by default), for the first set only.

check prints, for each build:
- whether the walkers share the same bias (zed and nker equal on all walkers at all times);
- the largest difference between the runs with and without the neighbor list (zero if exact,
  in serial walkers; round-off that grows along the trajectory with domain decomposition);
- transitions, and the final DeltaF from the STATE file of walker 0 and from reweighting all walkers.

env: GROMACS (default 2024.3), SETS (default 4), JOBS (cores, default nproc)
"""

import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from opesval.analysis import (  # noqa: E402
    deltaF,
    fes_from_state,
    force_shifted_kernels,
    make_grid,
    read_colvar,
    read_states,
)
from opesval.run import add_opes_keywords  # noqa: E402
from opesval.systems import SYSTEMS  # noqa: E402

OPT = ROOT / "opt"
MPI = OPT / "openmpi"
SYSTEM = SYSTEMS["ala2"]
METHODS = {"opes": "plumed-opes.dat", "explore": "plumed-explore.dat"}
LAYOUTS = {  # walkers, MPI ranks per walker, PLUMED threads, MD steps
    "walkers": (4, 1, 1, 2_500_000),
    "dd": (2, 2, 2, 500_000),
}


def plumed_version(label):
    out = subprocess.run([OPT / "plumed" / label / "bin" / "plumed", "info", "--version"], capture_output=True,
                         text=True, env=dict(os.environ, LD_LIBRARY_PATH=f"{OPT}/plumed/{label}/lib:{MPI}/lib")).stdout
    major, minor = out.strip().split(".")[:2]
    return int(major), int(minor)


def plumed_input(method, label, nolist):
    text = (SYSTEM.path / METHODS[method]).read_text()
    if nolist:
        text = re.sub(r"^\s*NLIST\s*\n", "", text, flags=re.M)
        if plumed_version(label) >= (2, 11):
            text = add_opes_keywords(text, "NLIST_OFF")
    return add_opes_keywords(text, "WALKERS_MPI")


def jobs(outdir, labels, sets):
    out = []
    for label in labels:
        for method in METHODS:
            for layout, (nwalk, _, _, _) in LAYOUTS.items():
                for s in range(sets if layout == "walkers" else 1):
                    for nolist in ([False, True] if s == 0 else [False]):
                        name = f"{label}/{method}/{layout}/{'nolist' if nolist else 'list'}-set{s}"
                        out.append((outdir / name, label, method, layout, s, nolist))
    return out


def run_one(job):
    d, label, method, layout, s, nolist = job
    nwalk, ranks, threads, nsteps = LAYOUTS[layout]
    if (d / "done").exists():
        return d, 0.0, "reused"
    d.mkdir(parents=True, exist_ok=True)
    text = plumed_input(method, label, nolist)
    dirs = []
    for w in range(nwalk):
        wd = d / f"w{w}"
        wd.mkdir(exist_ok=True)
        replica = nwalk * s + w
        (wd / "input.tpr").unlink(missing_ok=True)
        (wd / "input.tpr").symlink_to(SYSTEM.path / "tpr" / f"input.{replica}.tpr")
        (wd / "plumed.dat").write_text(text)
        dirs.append(wd.name)
    prefix = OPT / "plumed" / label
    env = dict(os.environ, PLUMED_KERNEL=f"{prefix}/lib/libplumedKernel.so", PLUMED_NUM_THREADS=str(threads),
               OMP_NUM_THREADS="1", PATH=f"{MPI}/bin:" + os.environ["PATH"],
               LD_LIBRARY_PATH=f"{MPI}/lib:" + os.environ.get("LD_LIBRARY_PATH", ""))
    gmx = OPT / "gromacs" / (os.environ.get("GROMACS", "2024.3") + "-mpi") / "bin" / "gmx_mpi"
    cmd = [MPI / "bin" / "mpirun", "-np", str(nwalk * ranks), "--oversubscribe", "--bind-to", "none",
           "--mca", "mpi_yield_when_idle", "1",  # do not spin while waiting, in case cores are shared
           "-x", "PLUMED_KERNEL", "-x", "PLUMED_NUM_THREADS", "-x", "OMP_NUM_THREADS", "-x", "LD_LIBRARY_PATH",
           gmx, "mdrun", "-multidir", *dirs, "-s", "input.tpr", "-plumed", "plumed.dat", "-deffnm", "md",
           "-nsteps", str(nsteps), "-ntomp", "1", "-pin", "off", "-nb", "cpu"]
    if ranks > 1:
        cmd += ["-dd", str(ranks), "1", "1"]
    t0 = time.time()
    with open(d / "out.log", "w") as f:
        proc = subprocess.run([str(c) for c in cmd], cwd=d, env=env, stdout=f, stderr=subprocess.STDOUT)
    wall = time.time() - t0
    if proc.returncode != 0:
        return d, wall, f"FAILED (exit {proc.returncode}), see {d / 'out.log'}"
    sha = re.search(r"Version: \S+ \(git: ([0-9a-f]+)\)", (d / "w0" / "md.log").read_text(errors="replace"))
    (d / "done").write_text(f"wall_time={wall:.1f}\nplumed_git={sha.group(1) if sha else '?'}\n")
    return d, wall, "done"


def run(outdir, labels):
    sets = int(os.environ.get("SETS", 4))
    todo = jobs(outdir, labels, sets)
    cores = int(os.environ.get("JOBS", os.cpu_count()))
    # schedule by cores: each job takes walkers*ranks*threads of them
    todo.sort(key=lambda j: -np.prod(LAYOUTS[j[3]][:3]))
    print(f"{len(todo)} jobs on {cores} cores")
    with ThreadPoolExecutor(max_workers=max(1, cores // 4)) as pool:
        for d, wall, status in pool.map(run_one, todo):
            print(f"{d.relative_to(outdir)}: {status} ({wall:.0f}s)", flush=True)


def walker_colvars(d):
    # with -multidir PLUMED adds the walker index as suffix: COLVAR.0, COLVAR.1, ...
    walkers = sorted(d.glob("w[0-9]*"), key=lambda w: int(w.name[1:]))
    return [read_colvar(w / f"COLVAR.{w.name[1:]}") for w in walkers]


def check(outdir, labels):
    a = SYSTEM.analysis
    grid = make_grid(a.grid)
    in_B = a.in_B(np.meshgrid(*grid, indexing="ij")[0])
    kbt = SYSTEM.kbt
    for label in labels:
        for method in METHODS:
            for layout in LAYOUTS:
                runs = sorted(d for d in (outdir / label / method / layout).glob("list-set*") if (d / "done").exists())
                if not runs:
                    continue
                shared, trans, dF_fes, dF_rew = True, [], [], []
                for d in runs:
                    cols = walker_colvars(d)
                    n = min(len(c) for c in cols)
                    for key in ("opes.zed", "opes.nker"):
                        ref = cols[0][key].to_numpy()[:n]
                        shared &= all(np.array_equal(c[key].to_numpy()[:n], ref) for c in cols[1:])
                    for c in cols:
                        phi = c["phi"].to_numpy()
                        state = np.where(a.core_A(phi), 0, np.where(a.core_B(phi), 1, -1))
                        s = state[state >= 0]
                        trans.append(np.count_nonzero(s[1:] != s[:-1]))
                    states = read_states(d / "w0" / "STATE")
                    fs = force_shifted_kernels(d / "w0")
                    dF_fes.append(deltaF(fes_from_state(states[-1], grid, fs), in_B) * kbt)
                    half = [c.iloc[len(c) // 2:] for c in cols]
                    phi = np.concatenate([c["phi"].to_numpy() for c in half])
                    lw = np.concatenate([c["opes.bias"].to_numpy() for c in half]) / kbt
                    inB = a.in_B(phi)
                    dF_rew.append((np.logaddexp.reduce(lw[~inB]) - np.logaddexp.reduce(lw[inB])) * kbt)
                line = (f"{label:14s} {method:8s} {layout:8s} sets={len(runs)} shared_bias={'yes' if shared else 'NO'}"
                        f" transitions/walker={np.mean(trans):.0f} dF_state={np.mean(dF_fes):.2f}"
                        f"+-{np.std(dF_fes) / np.sqrt(len(runs)):.2f} dF_rew(2nd half)={np.mean(dF_rew):.2f}"
                        f"+-{np.std(dF_rew) / np.sqrt(len(runs)):.2f}")
                d0, d1 = outdir / label / method / layout / "list-set0", outdir / label / method / layout / "nolist-set0"
                if (d0 / "done").exists() and (d1 / "done").exists():
                    diffs, same = [], True
                    for c0, c1 in zip(walker_colvars(d0), walker_colvars(d1)):
                        n = min(len(c0), len(c1))
                        cols = [k for k in ("phi", "psi", "opes.bias", "opes.zed", "opes.nker") if k in c1]
                        x0, x1 = c0[cols].to_numpy()[:n], c1[cols].to_numpy()[:n]
                        same &= np.array_equal(x0, x1)
                        first = np.flatnonzero(np.any(x0 != x1, axis=1))
                        diffs.append(c0["time"].iloc[first[0]] if len(first) else np.inf)
                    line += f" | list vs nolist: {'identical' if same else f'first difference at t={min(diffs):.0f} ps'}"
                print(line)


if __name__ == "__main__":
    if len(sys.argv) < 4 or sys.argv[1] not in ("run", "check"):
        sys.exit(__doc__)
    outdir = Path(sys.argv[2]).resolve()
    (run if sys.argv[1] == "run" else check)(outdir, sys.argv[3:])
