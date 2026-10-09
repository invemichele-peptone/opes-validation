"""Figures for docs/kernel-cutoff.md.

usage: uv run python docs/make_figures.py
Figures 1 and 2 are analytic. Figure 3 needs the `transient` runs of ala2 for the
PLUMED labels v2.9-base, pr1463 and alt-fshift (see the document).
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "docs" / "figures"

# one color per version, in every figure (validated categorical slots 1-3)
COLOR = {"old": "#2a78d6", "PR": "#eb6834", "alt": "#1baf7a"}
NAME = {"old": "old (v2.9)", "PR": "PR #1463", "alt": "alternative"}
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


if __name__ == "__main__":
    figure1()
    figure2()
    figure3()
