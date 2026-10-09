"""Figures for docs/kernel-cutoff.md and docs/kernel-merging.md.

usage: uv run python docs/make_figures.py [figure ...]
The kernel-shape figures are analytic, the others need the runs described in the documents.
"""

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "docs" / "figures"

# one color per version, in every figure (validated categorical slots 1-3)
COLOR = {"old": "#2a78d6", "PR": "#eb6834", "alt": "#1baf7a", "PRsplit": "#eda100"}
NAME = {"old": "old (v2.9)", "PR": "PR #1463", "alt": "alternative", "PRsplit": "PR + SPLIT"}
RUN_SET = {"old": "v2.9-base", "PR": "pr1463", "alt": "alt-fshift",
           "PRsplit": "schemes+cutoff_scheme-pr+merge_scheme-split"}
INK, INK2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#8b8a85", "#e4e3df", "#fcfcfb"


def style(ax):
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


def label(ax, x, y, text, color, **kw):
    ax.annotate(text, (x, y), color=INK2, fontsize=8, **kw)
    ax.plot([], [], color=color)  # keep the legend entry color-only


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(OUT / name, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"written {OUT / name}")


# kernel shapes, as a function of the scaled distance r (units of sigma)
def shapes(rc):
    c = np.exp(-0.5 * rc**2)
    inside = lambda r, v: np.where(r < rc, v, 0.0)  # noqa: E731
    gauss = lambda r: np.exp(-0.5 * r**2)  # noqa: E731
    g = lambda r: inside(r, gauss(r) - c)  # energy-shifted value
    f = lambda r: inside(r, gauss(r) - c * (1 + 0.5 * (rc**2 - r**2)))  # force-shifted value
    dg = lambda r: inside(r, r * gauss(r))  # -dg/dr
    df = lambda r: inside(r, r * (gauss(r) - c))  # -df/dr
    return c, gauss, g, f, dg, df


def figure1():
    """Kernel value and slope for the two truncation schemes, with an exaggerated cutoff."""
    rc = 2.5
    c, gauss, g, f, dg, df = shapes(rc)
    r = np.linspace(0, 3.2, 2000)
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8), facecolor=SURFACE)

    ax = axs[0]
    ax.plot(r, gauss(r), color=MUTED, ls="--", lw=1.2, label="Gaussian (not truncated)")
    ax.plot(r, g(r), color=COLOR["PR"], lw=2, label="energy-shifted  g = Gaussian − c")
    ax.plot(r, f(r), color=COLOR["alt"], lw=2, label="force-shifted  f = g − c·(r$_c^2$ − r$^2$)/2")
    ax.set_xlim(1.0, 3.2)
    ax.set_ylim(-0.01, 0.62)
    ax.axvline(rc, color=INK2, lw=0.8, ls=":")
    ax.annotate("cutoff $r_c$", (rc, 0.3), xytext=(4, 0), textcoords="offset points", color=INK2, fontsize=8)
    ax.set_xlabel("distance from the kernel centre r (units of σ)")
    ax.set_ylabel("kernel value K(r)")
    ax.set_title("value", fontsize=10)

    ax = axs[1]
    ax.plot(r, r * gauss(r), color=MUTED, ls="--", lw=1.2, label="Gaussian (not truncated)")
    ax.plot(r, dg(r), color=COLOR["PR"], lw=2, label="slope of g")
    ax.plot(r, df(r), color=COLOR["alt"], lw=2, label="slope of f")
    ax.set_xlim(1.0, 3.2)
    ax.set_ylim(-0.01, 0.65)
    ax.axvline(rc, color=INK2, lw=0.8, ls=":")
    ax.annotate("cutoff $r_c$", (rc, 0.3), xytext=(4, 0), textcoords="offset points", color=INK2, fontsize=8)
    ax.annotate("g: jumps to 0", (rc, dg(rc - 1e-9)), xytext=(6, 4), textcoords="offset points", color=INK2,
                fontsize=8)
    ax.annotate("f: smoothly to 0", (rc, 0), xytext=(6, 6), textcoords="offset points", color=INK2, fontsize=8)
    ax.set_xlabel("distance from the kernel centre r (units of σ)")
    ax.set_ylabel("−dK/dr")
    ax.set_title("slope (gives the force)", fontsize=10)
    for ax in axs:
        style(ax)
        ax.legend(fontsize=8, frameon=False, labelcolor=INK2, loc="upper right")
    fig.suptitle(f"Truncated kernels, cutoff exaggerated to $r_c$ = {rc}σ (c = {c:.3f}) for visibility",
                 color=INK, fontsize=11)
    save(fig, "kernel-shapes.png")


def figure2():
    """OPES_METAD force from a single kernel, with the ala2 parameters of the Explore paper."""
    kT, barrier = 0.0083144626 * 300, 30.0
    gamma = barrier / kT
    a = (1 - 1 / gamma) * kT
    eps = np.exp(-barrier / (1 - 1 / gamma) / kT)
    rc = np.sqrt(2 * barrier / (1 - 1 / gamma) / kT)  # default KERNEL_CUTOFF, then c = eps
    c, gauss, g, f, dg, df = shapes(rc)
    r = np.linspace(3.5, 5.4, 4000)
    force = {
        "old": lambda hz: a * hz * df(r) / (hz * g(r) + eps),
        "PR": lambda hz: a * hz * dg(r) / (hz * g(r) + eps),
        "alt": lambda hz: a * hz * df(r) / (hz * f(r) + eps),
    }
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8), facecolor=SURFACE)
    for ax, hz, ymax in zip(axs, (0.1, 1.0), (10, 13)):
        for key in ("old", "PR", "alt"):
            ax.plot(r, force[key](hz), color=COLOR[key], lw=2, label=NAME[key])
        ax.axvline(rc, color=INK2, lw=0.8, ls=":")
        ax.annotate("cutoff", (rc, ymax * 0.6), xytext=(4, 0), textcoords="offset points", color=INK2, fontsize=8)
        i = np.searchsorted(r, rc) - 1
        ax.annotate(f"← PR jumps by {force['PR'](hz)[i]:.1f}", (rc, force["PR"](hz)[i]), xytext=(4, -3),
                    textcoords="offset points", color=INK2, fontsize=8)
        ax.set_xlim(r[0], r[-1])
        ax.set_ylim(0, ymax)
        ax.set_xlabel("distance from the kernel centre r (units of σ)")
        ax.set_ylabel("bias force −dV/dr (kJ/mol per σ)")
        ax.set_title(f"kernel weight h/Z = {hz:g}", fontsize=10)
        style(ax)
        ax.legend(fontsize=8, frameon=False, labelcolor=INK2, loc="lower left")
    fig.suptitle(f"Force of one OPES_METAD kernel, ala2 setup: BARRIER=30 kJ/mol, $r_c$={rc:.2f}σ, c=ε={eps:.1e}",
                 color=INK, fontsize=11)
    save(fig, "kernel-force.png")


def figure3():
    """Early transient of ala2, 400 independent replicas per version."""
    runs = ROOT / "runs"
    labels = {"old": "v2.9-base", "PR": "pr1463", "alt": "alt-fshift"}
    reps = range(100, 500)
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.8), facecolor=SURFACE)
    for ax, method, title in zip(axs[:2], ("opes", "explore"),
                                 ("OPES_METAD", "OPES_METAD_EXPLORE (control: change ≈ 1e-5)")):
        for key, lab in labels.items():
            dirs = [runs / lab / "transient/ala2" / method / f"{i:03d}" for i in reps]
            if not all((d / "analysis.npz").exists() for d in dirs):
                print(f"missing runs for {lab}/{method}, skipping figure 3")
                plt.close(fig)
                return
            x = np.array([np.load(d / "analysis.npz")["deltaF_fes"] for d in dirs])
            t = np.load(dirs[0] / "analysis.npz")["time"]
            m, s = x.mean(0), x.std(0, ddof=1) / np.sqrt(len(x))
            ax.plot(t, m, color=COLOR[key], lw=2, marker="o", ms=4, label=NAME[key])
            ax.fill_between(t, m - s, m + s, color=COLOR[key], alpha=0.18, lw=0)
            if method == "opes":
                dy = {"old": 7, "alt": 0, "PR": -7}[key]  # spread the labels of close curves
                ax.annotate(NAME[key], (t[-1], m[-1]), xytext=(6, dy - 3), textcoords="offset points", color=INK2,
                            fontsize=8)
            if method == "opes" and key == "old":
                tf = {k: np.array([json.load(open(runs / lb / "transient/ala2/opes" / f"{i:03d}" / "analysis.json"))
                                   ["t_first"] for i in reps]) for k, lb in labels.items()}
        ax.axhline(9.05, color=INK2, ls="--", lw=1, label="reference ΔF")
        ax.set_xlim(t[0] - 30, t[-1] + (110 if method == "opes" else 30))
        ax.legend(fontsize=8, frameon=False, labelcolor=INK2, loc="upper right")
        ax.set_xlabel("time (ps)")
        ax.set_ylabel("ΔF from the bias (kJ/mol)")
        ax.set_title(title, fontsize=10)
        style(ax)
    ax = axs[2]
    for key, v in tf.items():
        v = np.sort(v)
        ax.step(np.r_[v[0], v], np.r_[0, np.arange(1, len(v) + 1) / len(v)], where="post", color=COLOR[key], lw=2,
                label=f"{NAME[key]} (median {np.median(v):.0f} ps)")
    ax.set_xscale("log")
    ax.set_xlim(5, 1000)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("time (ps)")
    ax.set_ylabel("fraction of replicas")
    ax.set_title("OPES_METAD: first transition to basin B", fontsize=10)
    style(ax)
    ax.legend(fontsize=8, frameon=False, labelcolor=INK2, loc="lower right")
    fig.suptitle("ala2, first 0.5 ns: mean ± s.e.m. over 400 independent replicas per version", color=INK, fontsize=11)
    save(fig, "transient-ala2.png")


def gauss(x, h, c, s):
    return h * np.exp(-0.5 * ((x - c) / s) ** 2)


def merged(h1, c1, s1, h2, c2, s2, split):
    """Merge two 1D kernels, as in OPESmetad::mergeKernels."""
    m1, m2 = (h1 * s1, h2 * s2) if split else (h1, h2)
    w1, w2 = m1 / (m1 + m2), m2 / (m1 + m2)
    c = c1 + w2 * (c2 - c1)
    s = np.sqrt(w1 * s1**2 + w2 * s2**2 + w1 * w2 * (c2 - c1) ** 2)
    return ((m1 + m2) / s if split else m1 + m2), c, s


def figure_merge_example():
    """Merging two kernels with the PLAIN and SPLIT rules."""
    x = np.linspace(-4, 4, 2000)
    cases = [((1.0, -0.5, 1.0), (1.0, 0.5, 1.0), "equal kernels 1σ apart"),
             ((1.0, -0.3, 1.0), (2.0, 0.3, 0.5), "old wide kernel + new narrow kernel")]
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8), facecolor=SURFACE)
    for ax, (k1, k2, title) in zip(axs, cases):
        exact = gauss(x, *k1) + gauss(x, *k2)
        mass = lambda h, c, s: h * s  # noqa: E731, up to sqrt(2 pi)
        ax.plot(x, gauss(x, *k1), color=MUTED, lw=1, ls=":")
        ax.plot(x, gauss(x, *k2), color=MUTED, lw=1, ls=":", label="the two kernels")
        ax.plot(x, exact, color=INK2, lw=1.5, ls="--", label=f"their sum (mass {mass(*k1) + mass(*k2):.2f})")
        for key, split in (("PR", False), ("PRsplit", True)):
            h, c, s = merged(*k1, *k2, split)
            name = "PLAIN: heights summed" if not split else "SPLIT: weights summed"
            ax.plot(x, gauss(x, h, c, s), color=COLOR[key], lw=2,
                    label=f"{name} (mass {mass(h, c, s):.2f}, centre {c:+.2f})")
        ax.set_xlim(x[0], x[-1])
        ax.set_ylim(0, 1.6 * exact.max())  # room for the legend
        ax.set_xlabel("CV")
        ax.set_ylabel("kernel value")
        ax.set_title(title, fontsize=10)
        style(ax)
        ax.legend(fontsize=7.5, frameon=False, labelcolor=INK2, loc="upper left")
    fig.suptitle("Merging two kernels (mass = height × σ, up to √(2π))", color=INK, fontsize=11)
    save(fig, "merge-example.png")


def mass_ratio(run, ncv):
    """KDE mass of the compressed kernels over that of all deposited kernels, at each stored STATE."""
    from opesval.analysis import read_states

    dep = np.loadtxt(run / "KERNELS", comments="#", ndmin=2)
    dep_mass = np.cumsum(dep[:, 1 + 2 * ncv] * np.prod(dep[:, 1 + ncv:1 + 2 * ncv], axis=1))
    out = []
    for st in read_states(run / "STATE"):
        k = st["kernels"]
        n = np.searchsorted(dep[:, 0], st["time"] + 1e-9)
        out.append((st["time"], (k[:, 1 + 2 * ncv] * np.prod(k[:, 1 + ncv:1 + 2 * ncv], axis=1)).sum() / dep_mass[n - 1]))
    return np.array(out)


def figure_mass_inflation():
    """How much the PLAIN merge inflates the mass of the KDE during the quick-tier runs."""
    runs = ROOT / "runs"
    cases = [("wolfe-quapp", "opes", 1), ("mueller", "opes", 1), ("mueller", "explore", 1),
             ("ala2", "opes_jpcl", 2), ("ala2", "opes", 2), ("ala2", "explore", 2)]
    fig, axs = plt.subplots(2, 3, figsize=(12, 6.2), facecolor=SURFACE)
    for ax, (system, method, ncv) in zip(axs.flat, cases):
        for key in ("PR", "PRsplit"):
            dirs = sorted((runs / RUN_SET[key] / "quick" / system / method).glob("[0-9][0-9][0-9]"))
            if not dirs:
                print(f"missing runs for {key} {system}/{method}, skipping")
                plt.close(fig)
                return
            ratios = [mass_ratio(d, ncv) for d in dirs]
            n = min(len(r) for r in ratios)
            t = ratios[0][:n, 0]
            r = np.array([x[:n, 1] for x in ratios])
            m, e = r.mean(0), r.std(0, ddof=1) / np.sqrt(len(r))
            ax.plot(t, m, color=COLOR[key], lw=2, label=NAME[key].replace("PR #1463", "PR (PLAIN merge)"))
            ax.fill_between(t, m - e, m + e, color=COLOR[key], alpha=0.2, lw=0)
            ax.annotate(f"{m[-1]:.2f}", (t[-1], m[-1]), xytext=(4, -3), textcoords="offset points", color=INK2,
                        fontsize=8)
        unit = " (ps)" if system == "ala2" else ""
        ax.set_xlim(0, t[-1] * 1.12)
        ax.set_ylim(0.9, None)
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
        ax.set_xlabel(f"time{unit}")
        ax.set_ylabel("mass compressed / deposited")
        ax.set_title(f"{system} / {method}", fontsize=10)
        style(ax)
    axs.flat[0].legend(fontsize=8, frameon=False, labelcolor=INK2, loc="lower right")
    fig.suptitle("KDE mass after compression, relative to the deposited kernels (mean ± s.e.m. over replicas)",
                 color=INK, fontsize=11)
    save(fig, "merge-mass.png")


def figure_ala4():
    """Exploration of ala4: number of the 8 states visited, quick tier."""
    runs = ROOT / "runs"
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8), facecolor=SURFACE)
    for ax, method, title in zip(axs, ("explore", "opes"), ("OPES_METAD_EXPLORE", "OPES_METAD")):
        for key in ("old", "PR", "PRsplit"):
            dirs = sorted((runs / RUN_SET[key] / "quick/ala4" / method).glob("[0-9][0-9][0-9]"))
            if not dirs:
                print(f"missing ala4 runs for {key}, skipping")
                plt.close(fig)
                return
            data = [np.load(d / "analysis.npz") for d in dirs]
            n = min(len(d["time"]) for d in data)
            t = data[0]["time"][:n]
            v = np.array([d["visited"][:n] for d in data])
            m, e = v.mean(0), v.std(0, ddof=1) / np.sqrt(len(v))
            ax.plot(t, m, color=COLOR[key], lw=2, label=NAME[key])
            ax.fill_between(t, m - e, m + e, color=COLOR[key], alpha=0.18, lw=0)
        ax.set_xscale("log")
        ax.set_xlim(10, t[-1])
        ax.set_ylim(0.8, 8.2)
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, _: f"{x:g}"))
        ax.set_xlabel("time (ps)")
        ax.set_ylabel("states visited (of 8)")
        ax.set_title(title, fontsize=10)
        style(ax)
        ax.legend(fontsize=8, frameon=False, labelcolor=INK2, loc="upper left")
    fig.suptitle("ala4, psi CVs: exploration of the 8 phi states (10 replicas, mean ± s.e.m.)", color=INK, fontsize=11)
    save(fig, "merge-ala4.png")


def figure_transient_split():
    """Early transient of ala2 with PR + SPLIT, 400 independent replicas per version."""
    runs = ROOT / "runs"
    keys = ("old", "PR", "PRsplit")
    reps = range(100, 500)
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.8), facecolor=SURFACE)
    tf = {}
    for ax, method in zip(axs[:2], ("opes", "explore")):
        for key in keys:
            dirs = [runs / RUN_SET[key] / "transient/ala2" / method / f"{i:03d}" for i in reps]
            if not all((d / "analysis.npz").exists() for d in dirs):
                print(f"missing transient runs for {key}, skipping")
                plt.close(fig)
                return
            x = np.array([np.load(d / "analysis.npz")["deltaF_fes"] for d in dirs])
            t = np.load(dirs[0] / "analysis.npz")["time"]
            m, e = x.mean(0), x.std(0, ddof=1) / np.sqrt(len(x))
            ax.plot(t, m, color=COLOR[key], lw=2, marker="o", ms=4, label=NAME[key])
            ax.fill_between(t, m - e, m + e, color=COLOR[key], alpha=0.18, lw=0)
            if method == "opes":
                ax.annotate(NAME[key], (t[-1], m[-1]), xytext=(6, {"old": 7, "PRsplit": 0, "PR": -7}[key] - 3),
                            textcoords="offset points", color=INK2, fontsize=8)
                tf[key] = np.array([json.load(open(d / "analysis.json"))["t_first"] for d in dirs])
        ax.axhline(9.05, color=INK2, ls="--", lw=1, label="reference ΔF")
        ax.set_xlim(t[0] - 30, t[-1] + (110 if method == "opes" else 30))
        ax.legend(fontsize=8, frameon=False, labelcolor=INK2, loc="upper right")
        ax.set_xlabel("time (ps)")
        ax.set_ylabel("ΔF from the bias (kJ/mol)")
        ax.set_title("OPES_METAD" if method == "opes" else "OPES_METAD_EXPLORE", fontsize=10)
        style(ax)
    ax = axs[2]
    for key, v in tf.items():
        v = np.sort(v)
        ax.step(np.r_[v[0], v], np.r_[0, np.arange(1, len(v) + 1) / len(v)], where="post", color=COLOR[key], lw=2,
                label=f"{NAME[key]} (median {np.median(v):.0f} ps)")
    ax.set_xscale("log")
    ax.set_xlim(5, 1000)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("time (ps)")
    ax.set_ylabel("fraction of replicas")
    ax.set_title("OPES_METAD: first transition to basin B", fontsize=10)
    style(ax)
    ax.legend(fontsize=8, frameon=False, labelcolor=INK2, loc="lower right")
    fig.suptitle("ala2, first 0.5 ns: mean ± s.e.m. over 400 independent replicas per version", color=INK, fontsize=11)
    save(fig, "merge-transient-ala2.png")


ZED_CASES = [("wolfe-quapp", "opes", 1), ("mueller", "opes", 1), ("mueller", "explore", 1),
             ("ala2", "opes_jpcl", 2), ("ala2", "opes", 2), ("ala2", "explore", 2),
             ("ala4", "opes", 3), ("ala4", "explore", 3)]


def load_zed(key, system, method, npoints=400):
    """Time and opes.zed of all quick-tier replicas, on a common subsampled time grid."""
    from opesval.analysis import read_colvar

    dirs = sorted((ROOT / "runs" / RUN_SET[key] / "quick" / system / method).glob("[0-9][0-9][0-9]"))
    if not dirs:
        return None, None
    cols = [read_colvar(d / "COLVAR")[["time", "opes.zed"]].to_numpy() for d in dirs]
    n = min(len(c) for c in cols)
    keep = np.unique(np.linspace(0, n - 1, min(n, npoints)).round().astype(int))
    return cols[0][keep, 0], np.array([c[keep, 1] for c in cols])


def figure_zed():
    """Evolution of the normalization Z in the quick-tier runs."""
    fig, axs = plt.subplots(2, 4, figsize=(15, 6.2), facecolor=SURFACE)
    for ax, (system, method, _) in zip(axs.flat, ZED_CASES):
        for key in ("old", "PR", "PRsplit"):
            t, z = load_zed(key, system, method)
            if t is None:
                print(f"missing runs for {key} {system}/{method}, skipping")
                plt.close(fig)
                return
            lz = np.log(z)  # Z is skewed across replicas: geometric mean, error from the s.e.m. of log Z
            m, e = lz.mean(0), lz.std(0, ddof=1) / np.sqrt(len(lz))
            ax.plot(t, np.exp(m), color=COLOR[key], lw=1.6, label=NAME[key])
            ax.fill_between(t, np.exp(m - e), np.exp(m + e), color=COLOR[key], alpha=0.18, lw=0)
        unit = " (ps)" if system.startswith("ala") else ""
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(t[1], t[-1])
        plain = matplotlib.ticker.FuncFormatter(lambda x, _: f"{x:g}")
        ax.xaxis.set_major_formatter(plain)
        ax.yaxis.set_major_locator(matplotlib.ticker.LogLocator(subs=(1.0, 2.0, 5.0)))
        ax.yaxis.set_major_formatter(plain)
        ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xlabel(f"time{unit}")
        ax.set_ylabel("Z")
        ax.set_title(f"{system} / {method}", fontsize=10)
        style(ax)
    axs.flat[0].legend(fontsize=8, frameon=False, labelcolor=INK2, loc="lower left")
    fig.suptitle("Normalization Z along the quick-tier runs (geometric mean ± s.e.m. over replicas)", color=INK,
                 fontsize=11)
    save(fig, "merge-zed.png")


def figure_zed_ratio():
    """Z(PR + SPLIT) / Z(PR) compared with the inverse of the mass inflation of the PLAIN merge."""
    fig, axs = plt.subplots(2, 3, figsize=(12, 6.2), facecolor=SURFACE)
    for ax, (system, method, ncv) in zip(axs.flat, ZED_CASES[:6]):
        t, zp = load_zed("PR", system, method)
        _, zs = load_zed("PRsplit", system, method)
        if t is None or zs is None:
            print(f"missing runs for {system}/{method}, skipping")
            plt.close(fig)
            return
        # replica by replica: paired runs share the adaptive sigma0, set before any merge, which fixes
        # the overall scale of Z (it varies a lot between replicas, much less between PR and PR + SPLIT)
        n = min(len(zp), len(zs))
        lr = np.log(zs[:n]) - np.log(zp[:n])
        d, e = lr.mean(0), lr.std(0, ddof=1) / np.sqrt(n)
        ratio = np.exp(d)
        ax.plot(t, ratio, color=COLOR["PRsplit"], lw=2, label="Z(PR + SPLIT) / Z(PR)")
        ax.fill_between(t, np.exp(d - e), np.exp(d + e), color=COLOR["PRsplit"], alpha=0.2, lw=0)
        dirs = sorted((ROOT / "runs" / RUN_SET["PR"] / "quick" / system / method).glob("[0-9][0-9][0-9]"))
        infl = [mass_ratio(d, ncv) for d in dirs]
        n = min(len(x) for x in infl)
        ti = infl[0][:n, 0]
        inv = np.exp(-np.log(np.array([x[:n, 1] for x in infl])).mean(0))
        ax.plot(ti, inv, color=INK2, lw=1.5, ls="--", label="1 / mass inflation of PR (PLAIN merge)")
        ax.axhline(1, color=MUTED, lw=0.8, ls=":")
        ax.annotate(f"{ratio[-1]:.2f}", (t[-1], ratio[-1]), xytext=(4, -10), textcoords="offset points",
                    color=INK2, fontsize=8)
        ax.annotate(f"{inv[-1]:.2f}", (ti[-1], inv[-1]), xytext=(4, 4), textcoords="offset points",
                    color=INK2, fontsize=8)
        unit = " (ps)" if system.startswith("ala") else ""
        ax.set_xlim(0, t[-1] * 1.12)
        ax.set_ylim(0.5, 1.2)
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
        ax.set_xlabel(f"time{unit}")
        ax.set_ylabel("ratio")
        ax.set_title(f"{system} / {method}", fontsize=10)
        style(ax)
    axs.flat[0].legend(fontsize=8, frameon=False, labelcolor=INK2, loc="lower left")
    fig.suptitle("Z with SPLIT relative to PLAIN, replica by replica, against the inverse of the mass inflation",
                 color=INK, fontsize=11)
    save(fig, "merge-zed-ratio.png")


# ala4 OPES_METAD_EXPLORE variants of the "many" tier: name -> (run set, color, merge rule, NO_ZED)
MANY = {
    "PR": ("schemes+cutoff_scheme-pr", COLOR["PR"], "PLAIN", False),
    "PR + MASS": ("schemes2+cutoff_scheme-pr+merge_scheme-mass", "#e87ba4", "MASS", False),
    "PR + SPLIT": ("schemes+cutoff_scheme-pr+merge_scheme-split", COLOR["PRsplit"], "SPLIT", False),
    "PR + SHAPE": ("schemes2+cutoff_scheme-pr+merge_scheme-shape", "#008300", "SHAPE", False),
    "PR, NO_ZED": ("schemes+cutoff_scheme-pr+no_zed", COLOR["PR"], "PLAIN", True),
    "PR + MASS, NO_ZED": ("schemes2+cutoff_scheme-pr+merge_scheme-mass+no_zed", "#e87ba4", "MASS", True),
    "PR + SPLIT, NO_ZED": ("schemes+cutoff_scheme-pr+merge_scheme-split+no_zed", COLOR["PRsplit"], "SPLIT", True),
    "PR + SHAPE, NO_ZED": ("schemes2+cutoff_scheme-pr+merge_scheme-shape+no_zed", "#008300", "SHAPE", True),
}


def load_many(run_set):
    dirs = sorted((ROOT / "runs" / run_set / "many/ala4/explore").glob("[0-9][0-9][0-9]"))
    dirs = [d for d in dirs if (d / "analysis.json").exists()]
    if len(dirs) < 50:
        return None
    vis = [np.load(d / "analysis.npz") for d in dirs]
    n = min(len(v["time"]) for v in vis)
    return {"time": vis[0]["time"][:n], "visited": np.array([v["visited"][:n] for v in vis]),
            "n_trans": np.array([json.load(open(d / "analysis.json"))["n_trans"] for d in dirs])}


def own_forces(name, t=5000.0, nrep=15):
    """Mean bias force at the visited points and number of kernels at time t, replaying the merge rule
    of a variant on its own trajectories (the first nrep replicas). Cached in reports/."""
    cache = ROOT / "reports" / "replay-ala4-own.json"
    data = json.load(open(cache)) if cache.exists() else {}
    if name in data:
        return np.array(data[name])
    from opesval.analysis import read_colvar
    from opesval.replay import bias, replay

    run_set, _, rule, no_zed = MANY[name]
    kT, periods, cvs = 0.0083144626181532 * 300, [2 * np.pi] * 3, ["psi1", "psi2", "psi3"]
    vals = []
    for run in sorted((ROOT / "runs" / run_set / "many/ala4/explore").glob("[0-9][0-9][0-9]"))[:nrep]:
        hd = {}
        for line in open(run / "KERNELS"):
            if not line.startswith("#!"):
                break
            key = line.split()[2] if line.startswith("#! SET") else None
            if key in ("biasfactor", "epsilon", "kernel_cutoff"):
                hd[key] = float(line.split()[3])
        a = (hd["biasfactor"] - 1) * kT
        dep = np.loadtxt(run / "KERNELS", comments="#", ndmin=2)
        ndep = np.searchsorted(dep[:, 0], t + 1e-9)
        col = read_colvar(run / "COLVAR")
        pts = col[col["time"] <= t][cvs].to_numpy()
        pts = pts[np.linspace(0, len(pts) - 1, min(800, len(pts))).astype(int)]
        k = replay(dep, 3, periods, rule, until=t)
        f = lambda p: bias(p, k, periods, hd["kernel_cutoff"], 1 + ndep, hd["epsilon"], a, no_zed)[0]  # noqa: E731
        g = np.zeros_like(pts)
        for i in range(3):
            dp, dm = pts.copy(), pts.copy()
            dp[:, i] += 1e-4
            dm[:, i] -= 1e-4
            g[:, i] = (f(dp) - f(dm)) / 2e-4
        vals.append((float(np.mean(np.linalg.norm(g, axis=1)) / kT), len(k[2])))
    data[name] = vals
    cache.parent.mkdir(exist_ok=True)
    json.dump(data, open(cache, "w"))
    return np.array(vals)


def figure_many():
    """ala4 explore with 50 replicas: exploration of each variant, against bias force and kernel number."""
    data = {name: load_many(v[0]) for name, v in MANY.items()}
    data = {k: v for k, v in data.items() if v is not None}
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.4), facecolor=SURFACE)
    ax = axs[0]
    for name, d in data.items():
        ax.plot(d["time"], d["visited"].mean(0), color=MANY[name][1], lw=1.8, ls="--" if MANY[name][3] else "-",
                label=name)
    ax.set_xscale("log")
    ax.set_xlim(10, d["time"][-1])
    ax.set_ylim(0.8, 8.2)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, _: f"{x:g}"))
    ax.set_xlabel("time (ps)")
    ax.set_ylabel("states visited (of 8), mean of 50")
    ax.set_title("exploration of the 8 φ states", fontsize=10)
    style(ax)
    ax.legend(fontsize=7.5, frameon=False, labelcolor=INK2, loc="upper left")
    for ax, col, xlabel in ((axs[1], 0, "bias force at visited points (kBT/rad)"),
                            (axs[2], 1, "compressed kernels")):
        for name, d in data.items():
            v = own_forces(name)[:, col]
            y, ye = d["n_trans"].mean(), d["n_trans"].std(ddof=1) / np.sqrt(len(d["n_trans"]))
            ax.errorbar(v.mean(), y, xerr=v.std(ddof=1) / np.sqrt(len(v)), yerr=ye, color=MANY[name][1],
                        marker="s" if MANY[name][3] else "o", ms=7, lw=1.5, capsize=2, ls="none")
            offset = (6, -13) if name == "PR, NO_ZED" else (6, 5)  # avoid the label of PR + MASS, NO_ZED
            ax.annotate(name, (v.mean(), y), xytext=offset, textcoords="offset points", color=INK2, fontsize=7.5)
        ax.set_xlabel(f"{xlabel}, replay of own runs at 5 ns")
        ax.set_ylabel("transitions in 10 ns, mean of 50")
        ax.margins(x=0.25)
        style(ax)
    axs[1].set_title("transitions vs bias force", fontsize=10)
    axs[2].set_title("transitions vs number of kernels", fontsize=10)
    fig.suptitle("ala4, OPES_METAD_EXPLORE with ψ CVs, 50 replicas × 10 ns per variant (dashed / squares: NO_ZED)",
                 color=INK, fontsize=11)
    save(fig, "merge-many-ala4.png")


UNCOMPRESSED = {  # name -> (run set, color); the reference first
    "uncompressed": ("schemes+cutoff_scheme-pr+compression_threshold-0", "#4a3aa7"),
    "PLAIN (PR)": ("pr1463", COLOR["PR"]),
    "MASS": ("schemes2+cutoff_scheme-pr+merge_scheme-mass", "#e87ba4"),
    "SPLIT": ("schemes+cutoff_scheme-pr+merge_scheme-split", COLOR["PRsplit"]),
}
UNCOMPRESSED_MANY = {"uncompressed": "schemes+cutoff_scheme-pr+compression_threshold-0", "PLAIN (PR)": "schemes+cutoff_scheme-pr",
                     "MASS": "schemes2+cutoff_scheme-pr+merge_scheme-mass",
                     "SPLIT": "schemes+cutoff_scheme-pr+merge_scheme-split"}


def load_scalars(run_set, tier, system, method):
    out = {}
    for d in sorted((ROOT / "runs" / run_set / tier / system / method).glob("[0-9][0-9][0-9]")):
        try:
            if json.load(open(d / "run.json"))["status"] == "done":
                out[int(d.name)] = json.load(open(d / "analysis.json"))
        except FileNotFoundError:
            pass
    return out


def figure_uncompressed():
    """Every merge rule against the uncompressed reference: transitions and free energy error."""
    cases = [("wolfe-quapp", "opes"), ("mueller", "opes"), ("mueller", "explore"), ("ala2", "opes_jpcl"),
             ("ala2", "opes"), ("ala2", "explore"), ("ala4", "opes"), ("ala4", "explore")]
    rows = []  # (label, {name: scalars by replica})
    for system, method in cases:
        data = {name: load_scalars(v[0], "quick", system, method) for name, v in UNCOMPRESSED.items()}
        if not data["uncompressed"]:
            continue
        common = set.intersection(*[set(d) for d in data.values()])
        rows.append((f"{system} / {method} ({len(common)})", {n: [d[i] for i in sorted(common)] for n, d in data.items()},
                     system))
    data = {n: list(load_scalars(rs, "many", "ala4", "explore").values()) for n, rs in UNCOMPRESSED_MANY.items()}
    if all(len(d) == 50 for d in data.values()):
        rows.append(("ala4 / explore, 10 ns (50)", data, "ala4"))
    if not rows:
        print("no uncompressed runs, skipping")
        return
    kT = {"ala2": 0.0083144626181532 * 300}
    fig, axs = plt.subplots(1, 2, figsize=(13, 0.55 * len(rows) + 1.8), facecolor=SURFACE)
    names = list(UNCOMPRESSED)
    offsets = dict(zip(names, np.linspace(-0.27, 0.27, len(names))))
    for y, (label, data, system) in enumerate(rows[::-1]):
        ref = np.array([r["n_trans"] for r in data["uncompressed"]], float)
        for name in names[1:]:
            x = np.array([r["n_trans"] for r in data[name]], float)
            ratio = x.mean() / ref.mean()
            err = ratio * np.sqrt((x.std(ddof=1) / np.sqrt(len(x)) / x.mean()) ** 2
                                  + (ref.std(ddof=1) / np.sqrt(len(ref)) / ref.mean()) ** 2)
            axs[0].errorbar(ratio, y + offsets[name], xerr=err, color=UNCOMPRESSED[name][1], marker="o", ms=6,
                            capsize=2, ls="none", label=name if y == 0 else None)
        if "rmse_deltaF_fes" in data["uncompressed"][0]:
            for name in names:
                x = np.array([r["rmse_deltaF_fes"] for r in data[name]], float) / kT.get(system, 1.0)
                axs[1].errorbar(x.mean(), y + offsets[name], xerr=x.std(ddof=1) / np.sqrt(len(x)),
                                color=UNCOMPRESSED[name][1], marker="o", ms=6, capsize=2, ls="none",
                                label=name if y == len(rows) - 1 or (y == 1 and "rmse" not in str(rows)) else None)
    labels = [r[0] for r in rows[::-1]]
    for ax in axs:
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels(labels, fontsize=8, color=INK2)
        ax.set_ylim(-0.6, len(rows) - 0.4)
        style(ax)
    axs[0].axvline(1, color=MUTED, lw=1, ls=":")
    axs[0].set_xlabel("transitions relative to the uncompressed runs")
    axs[0].set_title("exploration", fontsize=10)
    axs[0].legend(fontsize=8, frameon=False, labelcolor=INK2, loc="lower right")
    axs[1].set_xscale("log")
    axs[1].xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, _: f"{x:g}"))
    axs[1].set_xlabel("RMS error of ΔF in the 2nd half (kBT, log scale)")
    axs[1].set_title("free energy accuracy (two-basin systems)", fontsize=10)
    handles = [matplotlib.lines.Line2D([], [], color=UNCOMPRESSED[n][1], marker="o", ls="none") for n in names]
    axs[1].legend(handles, names, fontsize=8, frameon=False, labelcolor=INK2, loc="lower right")
    fig.suptitle("Merge rules against the uncompressed method (COMPRESSION_THRESHOLD=0); quick tier, replicas in brackets",
                 color=INK, fontsize=11)
    save(fig, "merge-uncompressed.png")


FIGURES = {f.__name__: f for f in (figure1, figure2, figure3, figure_merge_example, figure_mass_inflation,
                                   figure_ala4, figure_transient_split, figure_zed, figure_zed_ratio, figure_many,
                                   figure_uncompressed)}

if __name__ == "__main__":
    for name in sys.argv[1:] or FIGURES:
        FIGURES[name]()
