# OPES validation

Check that changes to PLUMED do not alter the qualitative behaviour of `OPES_METAD` and
`OPES_METAD_EXPLORE`, by rerunning the test systems of the original papers with two PLUMED
versions and comparing the results statistically.

- OPES paper: J. Phys. Chem. Lett. 11, 2731 (2020), <https://github.com/invemichele/opes>
- Explore paper: J. Chem. Theory Comput. 18, 3988 (2022), <https://github.com/invemichele/OPES-explore>

Systems: a modified Wolfe-Quapp and the Mueller-Brown 2D potentials (Langevin dynamics with
`plumed ves_md_linearexpansion`), and alanine di- and tetra-peptide in vacuum (GROMACS).
See [systems/README.md](systems/README.md) for the provenance of every input and the differences
from the papers.

## Quick start

Requirements: a C++ compiler, `make`, `cmake`, `git`, `curl`, and Python ≥ 3.10.

```bash
# python environment for the analysis
uv venv && uv pip install -r requirements.txt && source .venv/bin/activate
# or: python3 -m venv .venv && .venv/bin/pip install -r requirements.txt

# PLUMED versions to compare: build/build_plumed.sh <label> <ref> [repo]
build/build_plumed.sh v2.9-base d650dbd5a52ea8954950f823da5214f54bf837b5
build/build_plumed.sh pr1463 pr/1463

# GROMACS, patched in runtime mode with the oldest PLUMED to be tested
build/build_gromacs.sh 2024.3 v2.9-base

# run both versions, interleaved, on all cores
python -m opesval run -p v2.9-base pr1463 -t quick
# report in reports/v2.9-base_vs_pr1463/quick/report.md
python -m opesval compare v2.9-base pr1463 -t quick
```

Everything is installed in `opt/`, runs go in `runs/`, reports in `reports/` (all git-ignored).
Completed runs are reused, so `run` can be interrupted and restarted, or extended to more systems.

## How it works

- **One GROMACS for every PLUMED.** GROMACS is patched with `plumed patch --runtime`, so the
  PLUMED kernel is loaded at run time from `PLUMED_KERNEL`. The runner sets it for each job and
  checks the git hash printed in the log. Patch with the oldest PLUMED to be compared: newer
  kernels remain compatible with older patches. GROMACS 2024.3 has a patch in PLUMED v2.9, v2.10
  and master.
- **Identical inputs.** Replica *i* uses the same starting configuration and seeds for every
  PLUMED version. The trajectories still decorrelate after a short time, so runs are compared
  as distributions over replicas, not one by one.
- **Observables.**
  - ΔF between the basins, from the bias (STATE) and from reweighting, as a function of time.
  - Time to the first transition and number of transitions.
  - For ala4, the number of metastable states visited over time.
  - Final `rct`, `zed`, `neff` and number of kernels.
- **Statistics.** Mean ± s.e.m. over replicas, Welch's t-test and Mann-Whitney U. The report
  marks metrics with p < 0.01. With ~100 tests, about one is expected to be marked by chance;
  look for consistent patterns, not isolated values.

## Commands

```
python -m opesval list                          # systems, methods, tiers
python -m opesval run -p LABEL [LABEL ...]      # -s systems -m methods -t tier -r replicas -j jobs
python -m opesval analyze -p LABEL              # redo the analysis of completed runs
python -m opesval compare A B -t TIER           # write the report (--strict: exit 1 if anything is marked)
```

Use `--set NAME` with `run` to keep runs of the same PLUMED apart, e.g. with a different GROMACS
(`--gromacs 2025.0 --set v2.9-base_gmx2025`). Optional methods (6D ala4 from the OPES paper) run
only with `--optional` or when selected with `-m`.

## Tiers and cost

| tier | wolfe-quapp | mueller | ala2 | ala4 | wall time, 2 PLUMED versions on 20 cores |
|---|---|---|---|---|---|
| `smoke` | 4 × 4·10⁵ | 4 × 10⁶ | 4 × 2.5·10⁵ | 4 × 2.5·10⁵ | < 1 min |
| `quick` | 100 × 4·10⁶ | 25 × 2·10⁷ | 10 × 5·10⁶ | 10 × 10⁷ | ~45 min (estimate) |
| `full` (as in the papers) | 100 × 2·10⁷ | 25 × 2·10⁸ | 10 × 5·10⁶ | 10 × 10⁷ | ~4-5 h (estimate) |

Replicas × MD steps, per method; mueller dominates the cost. `smoke` only checks that everything
runs. It has no statistical power.

## Regenerating the GROMACS inputs

The `.tpr` files in `systems/*/tpr/` were generated with GROMACS 2024.3 and can be read by any
newer version. To regenerate them, e.g. after changing `md.mdp`:

```bash
build/make_tprs.sh ala2 2024.3
build/make_tprs.sh ala4 2024.3
```

## Automation

`.github/workflows/validate.yml` runs the same steps on GitHub Actions (manual trigger, with the
two PLUMED refs as inputs) and uploads the report as an artifact. GitHub-hosted runners have
4 cores, so only the `smoke` tier is practical there; `quick` and `full` need a self-hosted runner.
