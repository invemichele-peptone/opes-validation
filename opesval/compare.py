"""Statistical comparison of two run sets: A (reference) and B (candidate).

The trajectories of A and B decorrelate after a few thousand steps whatever the
change, so the comparison is between distributions over independent replicas.
"""

import json
import time
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker as ticker  # noqa: E402
import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

from .analysis import load_analysis  # noqa: E402
from .run import RUNS  # noqa: E402
from .systems import TwoBasins  # noqa: E402

P_FLAG = 0.01  # flag differences with p-value below this
CHECKPOINTS = (0.1, 0.25, 0.5, 1.0)  # fractions of the run where DeltaF(t) is tested

# validated categorical slots 1-2 of the dataviz reference palette, plus ink/grid tokens
COLOR = {"A": "#2a78d6", "B": "#eb6834"}
INK, INK2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#8b8a85", "#e4e3df", "#fcfcfb"


def collect(run_set, tier, system, method):
    """Completed replicas of a run set: {replica: (series, scalars, run_info)}."""
    out = {}
    for d in sorted((RUNS / run_set / tier / system.name / method.name).glob("[0-9][0-9][0-9]")):
        try:
            with open(d / "run.json") as f:
                info = json.load(f)
            if info["status"] == "done":
                out[int(d.name)] = (*load_analysis(d), info)
        except FileNotFoundError:
            continue
    return out


def stack(runs, kind):
    n = min(len(r[0][kind]) for r in runs)
    return runs[0][0]["time"][:n], np.array([r[0][kind][:n] for r in runs])


def mean_sem(x):
    x = np.asarray(x, dtype=float)
    n = np.sum(np.isfinite(x), axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.nanmean(x, axis=0), np.nanstd(x, axis=0, ddof=1) / np.sqrt(n)


def welch(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2 or (np.ptp(a) == 0 and np.ptp(b) == 0):
        return np.nan
    with warnings.catch_warnings():  # nearly identical samples, e.g. all replicas visited all states
        warnings.simplefilter("ignore", RuntimeWarning)
        return float(stats.ttest_ind(a, b, equal_var=False).pvalue)


def mwu(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2 or np.all(np.concatenate([a, b]) == a[0]):
        return np.nan
    return float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)


def checkpoints(n):
    """Indices of the time series where DeltaF(t) is tested."""
    return sorted({min(n - 1, max(0, int(round(frac * n)) - 1)) for frac in CHECKPOINTS})


def censored(values, t_end):
    """Times to an event that may not have happened (inf): rank them after every observed time."""
    v = np.asarray(values, float)
    return np.where(np.isfinite(v), v, 2 * t_end)


def fmt(x, digits=3):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "-"
    return f"{x:.{digits}g}"


def compare_method(system, method, runs_a, runs_b):
    """Returns the metric rows and the data needed for plotting."""
    common = sorted(set(runs_a) & set(runs_b))
    ra, rb = [runs_a[i] for i in common], [runs_b[i] for i in common]
    rows = []
    plot = {"n": len(common)}

    def add(metric, a, b, kind="mean"):
        if kind == "time":
            t_end = max(r[1]["t_end"] for r in ra + rb)
            ca, cb = censored(a, t_end), censored(b, t_end)
            med = lambda v: fmt(np.median(v)) if np.median(v) < t_end else "not reached"  # noqa: E731
            row = {"metric": metric, "A": f"median {med(ca)} ({np.isfinite(a).sum()}/{len(a)} reached)",
                   "B": f"median {med(cb)} ({np.isfinite(b).sum()}/{len(b)} reached)",
                   "p_welch": np.nan, "p_mwu": mwu(ca, cb)}
        else:
            (ma, sa), (mb, sb) = mean_sem(a), mean_sem(b)
            row = {"metric": metric, "A": f"{fmt(ma)} ± {fmt(sa, 2)}", "B": f"{fmt(mb)} ± {fmt(sb, 2)}",
                   "p_welch": welch(a, b), "p_mwu": mwu(a, b)}
        ps = [p for p in (row["p_welch"], row["p_mwu"]) if np.isfinite(p)]
        row["p"] = min(ps) if ps else np.nan
        row["flag"] = bool(ps) and bool(row["p"] < P_FLAG)
        rows.append(row)

    get = lambda runs, key: np.array([r[1].get(key, np.nan) for r in runs], float)  # noqa: E731
    if isinstance(system.analysis, TwoBasins):
        plot["ref"] = system.analysis.ref_deltaF() * system.kbt
        for kind in ("deltaF_fes", "deltaF_rew"):
            ta, xa = stack(ra, kind)
            tb, xb = stack(rb, kind)
            n = min(len(ta), len(tb))
            plot[kind] = (ta[:n], xa[:, :n], xb[:, :n])
            for i in checkpoints(n):
                add(f"{kind}(t={fmt(ta[i])})", xa[:, i], xb[:, i])
            add(f"rmse_{kind} (2nd half)", get(ra, f"rmse_{kind}"), get(rb, f"rmse_{kind}"))
            # largest undershoot, e.g. a basin overfilled right after the first transition
            add(f"min {kind}", np.nanmin(xa, axis=1), np.nanmin(xb, axis=1))
        add("t_first (A to B)", get(ra, "t_first"), get(rb, "t_first"), kind="time")
        plot["t_first"] = (get(ra, "t_first"), get(rb, "t_first"))
    else:
        ta, xa = stack(ra, "visited")
        tb, xb = stack(rb, "visited")
        n = min(len(ta), len(tb))
        plot["visited"] = (ta[:n], xa[:, :n], xb[:, :n])
        for i in checkpoints(n):
            add(f"visited states (t={fmt(ta[i])})", xa[:, i], xb[:, i])
        add("t_all (8 states visited)", get(ra, "t_all"), get(rb, "t_all"), kind="time")
        plot["t_all"] = (get(ra, "t_all"), get(rb, "t_all"))
    add("n_trans", get(ra, "n_trans"), get(rb, "n_trans"))
    for key in ("final_nker", "final_zed", "final_rct", "final_neff"):
        add(key, get(ra, key), get(rb, key))
    plot["wall_time"] = (np.mean([r[2]["wall_time"] for r in ra]), np.mean([r[2]["wall_time"] for r in rb]))
    return rows, plot


def _style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.xaxis.label.set_color(INK2)
    ax.yaxis.label.set_color(INK2)
    ax.title.set_color(INK)
    if ax.get_xscale() == "log":
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        ax.xaxis.set_minor_formatter(ticker.NullFormatter())


def _band(ax, t, x, label, color):
    m, s = mean_sem(x)
    ax.plot(t, m, color=color, lw=1.5, label=label)
    ax.fill_between(t, m - s, m + s, color=color, alpha=0.2, lw=0)


def _cdf(ax, values, label, color, t_end):
    v = np.sort(np.asarray(values, float))
    v = v[np.isfinite(v)]
    if len(v):
        frac = np.arange(1, len(v) + 1) / len(values)
        ax.step(np.r_[v[0], v, t_end], np.r_[0, frac, frac[-1]], where="post", color=color, lw=1.5, label=label)


def _cdf_panel(ax, values_a, values_b, labels, t_first_sample, t_end, unit, title):
    _cdf(ax, values_a, labels["A"], COLOR["A"], t_end)
    _cdf(ax, values_b, labels["B"], COLOR["B"], t_end)
    reached = np.concatenate([values_a, values_b])
    reached = reached[np.isfinite(reached)]
    if not len(reached):
        ax.text(0.5, 0.5, "no replica reached it", transform=ax.transAxes, ha="center", color=INK2, fontsize=9)
    ax.set_xscale("log")
    ax.set_xlim(min([t_first_sample, *reached]), t_end)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel(f"time{unit}")
    ax.set_ylabel("fraction of replicas")
    ax.set_title(title, fontsize=10)


def plot_method(system, method, plot, labels, paper_dir, path):
    unit = f" ({system.time_unit})" if system.time_unit else ""
    eunit = "kJ/mol" if system.kbt != 1 else "kBT"
    if "ref" in plot:
        fig, axs = plt.subplots(1, 3, figsize=(12, 3.6), facecolor=SURFACE)
        ref = plot["ref"]
        for ax, kind, title in zip(axs, ("deltaF_fes", "deltaF_rew"), ("ΔF from the bias (STATE)", "ΔF from reweighting")):
            t, xa, xb = plot[kind]
            _band(ax, t, xa, labels["A"], COLOR["A"])
            _band(ax, t, xb, labels["B"], COLOR["B"])
            ax.axhline(ref, color=INK2, ls="--", lw=1, label="reference")
            paper_file = system.paper_data.get(method.name, {}).get(kind)
            if paper_file:
                d = np.loadtxt(paper_dir / paper_file)
                ax.plot(d[:, 0], d[:, 1], color=MUTED, ls=":", lw=1.5, label="paper (mean)")
            ax.set_ylim(ref - 8 * system.kbt, ref + 12 * system.kbt)
            ax.set_xlim(0, t[-1])
            ax.set_xlabel(f"time{unit}")
            ax.set_ylabel(f"ΔF ({eunit})")
            ax.set_title(title, fontsize=10)
        t = plot["deltaF_fes"][0]
        _cdf_panel(axs[2], *plot["t_first"], labels, t[0] / 100, t[-1], unit, "first transition A→B")
    else:
        fig, axs = plt.subplots(1, 2, figsize=(8, 3.6), facecolor=SURFACE)
        t, xa, xb = plot["visited"]
        ax = axs[0]
        _band(ax, t, xa, labels["A"], COLOR["A"])
        _band(ax, t, xb, labels["B"], COLOR["B"])
        ax.set_xscale("log")
        ax.set_xlim(t[1], t[-1])
        ax.set_ylim(0.8, 8.2)
        ax.set_xlabel(f"time{unit}")
        ax.set_ylabel("visited states")
        ax.set_title("states visited (of 8)", fontsize=10)
        _cdf_panel(axs[1], *plot["t_all"], labels, t[1], t[-1], unit, "all 8 states visited")
    for ax in axs:
        _style(ax)
    axs[0].legend(fontsize=8, frameon=False, labelcolor=INK2)
    fig.suptitle(f"{system.name} / {method.name}  ({plot['n']} replicas, mean ± s.e.m.)", color=INK, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=130, facecolor=SURFACE)
    plt.close(fig)


def _provenance(run_set, tier, pairs):
    plumed, gromacs = {}, set()
    for system, method in pairs:
        for *_, info in collect(run_set, tier, system, method).values():
            plumed = info["plumed"]
            if info.get("gromacs"):
                gromacs.add(info["gromacs"])
    return plumed, ", ".join(sorted(gromacs))


def compare(set_a, set_b, tier, pairs, outdir):
    outdir = Path(outdir)
    (outdir / "figures").mkdir(parents=True, exist_ok=True)
    labels = {"A": set_a, "B": set_b}
    summary = {"A": set_a, "B": set_b, "tier": tier, "date": time.strftime("%Y-%m-%d %H:%M"), "results": {}}
    sections, flagged, ntests = [], [], 0
    for system, method in pairs:
        runs_a = collect(set_a, tier, system, method)
        runs_b = collect(set_b, tier, system, method)
        if not runs_a or not runs_b:
            continue
        rows, plot = compare_method(system, method, runs_a, runs_b)
        key = f"{system.name}/{method.name}"
        fig = f"figures/{system.name}-{method.name}.png"
        plot_method(system, method, plot, labels, system.path / "paper", outdir / fig)
        summary["results"][key] = rows
        ntests += sum(np.isfinite(r["p_welch"]) or np.isfinite(r["p_mwu"]) for r in rows)
        flagged += [(key, r) for r in rows if r["flag"]]
        lines = [f"### {key}", "", f"{system.description}; method from the {method.paper} paper, "
                 f"`{method.plumed}`. {plot['n']} replicas in common, "
                 f"mean wall time {plot['wall_time'][0]:.0f} s (A) / {plot['wall_time'][1]:.0f} s (B).", ""]
        if "ref" in plot:
            lines += [f"Reference ΔF = {plot['ref']:.3g} {'kJ/mol' if system.kbt != 1 else 'kBT'}.", ""]
        lines += [f"![{key}]({fig})", "", "| metric | A | B | p Welch | p MWU | |", "|---|---|---|---|---|---|"]
        for r in rows:
            lines.append(f"| {r['metric']} | {r['A']} | {r['B']} | {fmt(r['p_welch'], 2)} | {fmt(r['p_mwu'], 2)} "
                         f"| {'**check**' if r['flag'] else ''} |")
        sections.append("\n".join(lines))

    if not sections:
        raise SystemExit(f"nothing to compare for tier '{tier}' between '{set_a}' and '{set_b}'")
    plumed_a, gmx_a = _provenance(set_a, tier, pairs)
    plumed_b, gmx_b = _provenance(set_b, tier, pairs)
    head = [
        f"# OPES validation: `{set_a}` (A) vs `{set_b}` (B)",
        "",
        f"Tier `{tier}`, generated {summary['date']}.",
        "",
        "| | A | B |",
        "|---|---|---|",
        f"| PLUMED ref | `{plumed_a.get('ref')}` | `{plumed_b.get('ref')}` |",
        f"| PLUMED commit | `{plumed_a.get('sha')}` | `{plumed_b.get('sha')}` |",
        f"| GROMACS | {gmx_a or '-'} | {gmx_b or '-'} |",
        "",
        "Each replica of A and B starts from the same inputs and random seeds, but trajectories "
        "decorrelate quickly, so only the distributions over replicas are compared. "
        "Values are mean ± standard error of the mean over replicas (medians for times to an event, "
        "where runs that never reach the event are ranked last). "
        "p-values: Welch's t-test and Mann-Whitney U, both two-sided. "
        f"Rows with p < {P_FLAG} are marked **check**; with {ntests} tests, about {ntests * P_FLAG:.1f} "
        "are expected to be marked by chance alone.",
        "",
        "## Summary",
        "",
    ]
    if flagged:
        head += [f"{len(flagged)} metrics marked:", ""]
        head += [f"- {k}: {r['metric']}: A {r['A']}, B {r['B']} (p = {fmt(r['p'], 2)})" for k, r in flagged]
    else:
        head += [f"No metric with p < {P_FLAG}."]
    head += ["", "## Results", ""]
    (outdir / "report.md").write_text("\n".join(head) + "\n" + "\n\n".join(sections) + "\n")
    with open(outdir / "summary.json", "w") as f:
        json.dump(summary, f, indent=1, default=float)
    print(f"report written to {outdir / 'report.md'}: {len(flagged)} of {ntests} tests with p < {P_FLAG}")
    return flagged
