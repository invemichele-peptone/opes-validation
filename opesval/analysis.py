"""Per-run analysis: free energy differences, transitions and exploration metrics."""

import json

import numpy as np
import pandas as pd

from .systems import Ala4Basins, TwoBasins

MAX_STATES = 200  # evaluate at most this many stored STATEs per run


def read_colvar(path):
    with open(path) as f:
        fields = f.readline().split()
    if fields[:2] != ["#!", "FIELDS"]:
        raise ValueError(f"no FIELDS header in {path}")
    data = pd.read_csv(path, sep=r"\s+", comment="#", header=None, names=fields[2:])
    return data.apply(pd.to_numeric, errors="coerce").dropna()


def read_states(path):
    """Read a STATE file written with STORE_STATES. Returns a list of dicts."""
    states = []
    with open(path) as f:
        for line in f:
            if line.startswith("#! FIELDS"):
                fields = line.split()[2:]
                state = {"fields": fields, "set": {}, "rows": []}
                states.append(state)
            elif line.startswith("#! SET"):
                _, _, key, value = line.split(maxsplit=3)
                state["set"][key] = value.strip()
            elif line.strip() and not line.startswith("#"):
                state["rows"].append(line)
    for state in states:
        state["kernels"] = np.loadtxt(state.pop("rows"), ndmin=2)
        state["time"] = state["kernels"][0, 0]
    return states


def fes_from_state(state, grid):
    """Free energy estimate (in kBT) of the bias encoded in a stored STATE.

    Same estimator as FES_from_State.py of the paper repositories.
    `grid` is a list of 1D arrays, one per CV; returns an array of shape [len(g) for g in grid].
    """
    s = state["set"]
    ncv = len(grid)
    kernels = state["kernels"]
    centers = kernels[:, 1 : 1 + ncv]
    sigmas = kernels[:, 1 + ncv : 1 + 2 * ncv]
    heights = kernels[:, 1 + 2 * ncv]
    explore = s["action"] == "OPES_METAD_EXPLORE_state"
    zed = float(s["zed"]) * float(s["counter"] if explore else s["sum_weights"])
    epsilon = float(s["epsilon"])
    val_at_cutoff = np.exp(-0.5 * float(s["kernel_cutoff"]) ** 2)
    scale = float(s["biasfactor"]) if explore else 1.0

    points = np.stack([g.ravel() for g in np.meshgrid(*grid, indexing="ij")], axis=1)
    prob = np.empty(len(points))
    chunk = max(1, 4_000_000 // max(1, len(heights)))
    for start in range(0, len(points), chunk):
        p = points[start : start + chunk, None, :]
        dist = p - centers[None, :, :]
        for i, g in enumerate(grid):
            period = getattr(g, "period", None)
            if period:
                dist[:, :, i] -= period * np.round(dist[:, :, i] / period)
        norm2 = np.sum((dist / sigmas[None, :, :]) ** 2, axis=2)
        prob[start : start + chunk] = np.sum(heights * np.maximum(np.exp(-0.5 * norm2) - val_at_cutoff, 0), axis=1)
    prob = prob / zed + epsilon
    return (-scale * np.log(prob)).reshape([len(g) for g in grid])


class PeriodicGrid(np.ndarray):
    """A 1D grid that knows its period."""

    def __new__(cls, values, period):
        obj = np.asarray(values).view(cls)
        obj.period = period
        return obj

    def __array_finalize__(self, obj):
        self.period = getattr(obj, "period", None)


def make_grid(spec):
    grid = []
    for lo, hi, nbins, periodic in spec:
        if periodic:
            grid.append(PeriodicGrid(np.linspace(lo, hi, nbins, endpoint=False), hi - lo))
        else:
            grid.append(np.linspace(lo, hi, nbins))
    return grid


def deltaF(fes, in_B):
    """F(B) - F(A), with fes in kBT and in_B a boolean mask of the same shape."""
    return np.logaddexp.reduce(-fes[~in_B]) - np.logaddexp.reduce(-fes[in_B])


def running_deltaF_rew(cv, bias_over_kbt, in_B, idx):
    """Reweighted F(B)-F(A) in kBT, using all samples up to each index in idx."""
    lw_A = np.logaddexp.accumulate(np.where(in_B(cv), -np.inf, bias_over_kbt))
    lw_B = np.logaddexp.accumulate(np.where(in_B(cv), bias_over_kbt, -np.inf))
    with np.errstate(invalid="ignore"):
        out = lw_A[idx] - lw_B[idx]
    out[~np.isfinite(out)] = np.nan
    return out


def transitions(time, cv, core_A, core_B):
    """First passage to the other basin and number of core-to-core transitions."""
    state = np.full(len(cv), -1)
    state[core_A(cv)] = 0
    state[core_B(cv)] = 1
    in_core = np.flatnonzero(state >= 0)
    if len(in_core) == 0:
        return {"start_basin": -1, "t_first": np.inf, "n_trans": 0}
    s = state[in_core]
    changes = np.flatnonzero(s[1:] != s[:-1]) + 1
    t_first = time[in_core[changes[0]]] if len(changes) else np.inf
    return {"start_basin": int(s[0]), "t_first": float(t_first), "n_trans": int(len(changes))}


def final_values(colvar):
    last = colvar.iloc[-1]
    return {
        f"final_{k[5:]}": float(last[k])
        for k in ("opes.rct", "opes.zed", "opes.neff", "opes.nker")
        if k in colvar.columns
    }


def analyze_two_basins(system, run_dir):
    a = system.analysis
    colvar = read_colvar(run_dir / "COLVAR")
    time = colvar["time"].to_numpy()
    cv = colvar[a.cv].to_numpy()
    states = read_states(run_dir / "STATE")
    if len(states) > MAX_STATES:
        states = [states[i] for i in np.linspace(0, len(states) - 1, MAX_STATES).round().astype(int)]
    t_state = np.array([s["time"] for s in states])

    grid = make_grid(a.grid)
    mesh = np.meshgrid(*grid, indexing="ij")
    in_B = a.in_B(mesh[0])  # basins are defined on the first CV
    dF_fes = np.array([deltaF(fes_from_state(s, grid), in_B) for s in states])

    idx = np.clip(np.searchsorted(time, t_state, side="right") - 1, 0, len(time) - 1)
    dF_rew = running_deltaF_rew(cv, colvar["opes.bias"].to_numpy() / system.kbt, a.in_B, idx)

    series = {
        "time": t_state,
        "deltaF_fes": dF_fes * system.kbt,
        "deltaF_rew": dF_rew * system.kbt,
    }
    scalars = transitions(time, cv, a.core_A, a.core_B)
    scalars.update(final_values(colvar))
    ref = a.ref_deltaF() * system.kbt
    half = t_state >= t_state[-1] / 2
    for kind in ("deltaF_fes", "deltaF_rew"):
        err = series[kind][half] - ref
        scalars[f"rmse_{kind}"] = float(np.sqrt(np.nanmean(err**2))) if np.any(np.isfinite(err)) else np.nan
        scalars[f"final_{kind}"] = float(series[kind][-1])
    scalars["t_end"] = float(time[-1])
    return series, scalars


def analyze_ala4(system, run_dir):
    a = system.analysis
    colvar = read_colvar(run_dir / "COLVAR")
    time = colvar["time"].to_numpy()
    phis = colvar[list(a.phis)].to_numpy()
    # hysteresis on each phi, as in get_ala4_basin.sh
    to_1 = (phis > a.to_1[0]) & (phis < a.to_1[1])
    to_0 = (phis > a.to_0[0]) & (phis < a.to_0[1])
    bits = np.zeros(phis.shape, dtype=int)
    current = np.zeros(phis.shape[1], dtype=int)
    for n in range(len(time)):
        current = np.where(current == 1, np.where(to_0[n], 0, 1), np.where(to_1[n], 1, 0))
        bits[n] = current
    basin = bits @ (1 << np.arange(phis.shape[1]))  # 0..7
    first_visit = np.full(8, np.inf)
    for b in range(8):
        hit = np.flatnonzero(basin == b)
        if len(hit):
            first_visit[b] = time[hit[0]]
    visited = (first_visit[None, :] <= time[:, None]).sum(axis=1)
    keep = np.unique(np.linspace(0, len(time) - 1, min(len(time), 2000)).round().astype(int))
    series = {"time": time[keep], "visited": visited[keep].astype(float)}
    scalars = {
        "t_all": float(first_visit.max()),
        "final_visited": int(visited[-1]),
        "n_trans": int(np.count_nonzero(basin[1:] != basin[:-1])),
        "t_end": float(time[-1]),
    }
    scalars.update(final_values(colvar))
    return series, scalars


def analyze_run(system, run_dir):
    if isinstance(system.analysis, TwoBasins):
        series, scalars = analyze_two_basins(system, run_dir)
    elif isinstance(system.analysis, Ala4Basins):
        series, scalars = analyze_ala4(system, run_dir)
    else:
        raise TypeError(f"no analysis for {system.name}")
    np.savez(run_dir / "analysis.npz", **series)
    with open(run_dir / "analysis.json", "w") as f:
        json.dump(scalars, f, indent=1)
    return series, scalars


def load_analysis(run_dir):
    with np.load(run_dir / "analysis.npz") as data:
        series = dict(data)
    with open(run_dir / "analysis.json") as f:
        scalars = json.load(f)
    return series, scalars
