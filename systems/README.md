# Systems

All inputs come from the repositories of the two papers, with the changes listed below.

- **OPES paper**: M. Invernizzi and M. Parrinello, *Rethinking Metadynamics: From Bias Potentials to Probability Distributions*, J. Phys. Chem. Lett. 11, 2731 (2020), <https://github.com/invemichele/opes>
- **Explore paper**: M. Invernizzi and M. Parrinello, *Exploration vs Convergence Speed in Adaptive-Bias Enhanced Sampling*, J. Chem. Theory Comput. 18, 3988 (2022), <https://github.com/invemichele/OPES-explore>
- **Masterclass**: PLUMED Masterclass 22.3, <https://github.com/invemichele/masterclass-22-03>

| system | paper | engine | methods | replicas × steps (full) |
|---|---|---|---|---|
| `wolfe-quapp` | OPES | `plumed ves_md_linearexpansion` | `opes` | 100 × 2·10⁷ |
| `mueller` | Explore | `plumed ves_md_linearexpansion` | `opes`, `explore` | 25 × 2·10⁸ |
| `ala2` | OPES (`opes_jpcl`), Explore (`opes`, `explore`) | GROMACS | `opes_jpcl`, `opes`, `explore` | 10 × 5·10⁶ |
| `ala4` | Explore (`opes`, `explore`), OPES (`opes_jpcl`, optional) | GROMACS | `opes`, `explore`, `opes_jpcl` | 10 × 10⁷ |

## Changes with respect to the papers

Common to all systems:
- STATE output (`STATE_WFILE`, `STORE_STATES`) is added where missing, to compute the running free energy estimate. It has no effect on the dynamics.
- Every replica is a separate single-core run with a fixed seed, instead of multi-replica MPI runs.

`wolfe-quapp` (OPES paper `model/`): seeds and initial positions of the 10×10 replicas are taken from `model/inputs/md_input.{0..9}.dat` (`replicas.dat`). `ref_fes.dat` is `model/ref_fes-model.data`.

`mueller` (Explore paper `mueller/`): seeds and initial positions from `mueller/inputs/input_md.{0..24}.dat` (`replicas.dat`). `ref_fes.grid` is `mueller/static/modelFES.grid`. Its ΔF agrees with direct numerical integration of the 2D potential within 2·10⁻⁴ kBT.

`ala2` and `ala4`: the papers ran in vacuum with `cutoff-scheme = group`, `pbc = no`, which GROMACS removed in 2020. Following the masterclass, the runs here use the Verlet scheme with a periodic box and 1 nm plain cut-offs. Free energies are therefore slightly different from the papers. `ref_fes2D.dat` is the masterclass reference for this setup.
- `gromacs/topol.top` and `gromacs/md.mdp` are taken from the masterclass (`inputs-ala2/`, `inputs-ala4/`). The mdp differs only in fixed seeds (`gen_seed = ld_seed = 1000 + replica`) and, for ala4, `nsteps` from the Explore paper.
- `gromacs/conf.{0..9}.gro` are the 10 starting configurations of the papers (`opes/ala2/inputs/input.*.tpr` and `OPES-explore/ala4/inputs/input.*.tpr`), placed in the masterclass box (3.02334 nm for ala2, 5 nm for ala4), e.g.:
  ```
  gmx editconf -f opes/ala2/inputs/input.0.tpr -o conf.0.gro -box 3.02334 -c
  ```
  All of them start in the phi < 0 basin(s).
- `tpr/` is generated with `build/make_tprs.sh <system>` using GROMACS 2024.3. Replica *i* starts from `conf.<i % 10>.gro` with seed 1000 + *i*. ala2 has 100 replicas, for the `early` tier. `gmx check` against the masterclass `input.tpr` shows only the expected differences: coordinates, velocities, seeds, and `rlist`, which GROMACS 2024 recomputes for the Verlet buffer.
- ala2: `STATE_WSTRIDE=500*100` for `opes_jpcl`, instead of 10000.
- ala4: `ala4/explore/plumed.dat` in the Explore paper repository reads `OPES_METAD`, apparently a copy of `ala4/opes/plumed.dat`. Here `plumed-explore.dat` uses `OPES_METAD_EXPLORE`. Only the "bad CVs, standard parameters" setup (`psi-std`) is included.

## Observables

- `deltaF_fes`: ΔF = F(B) − F(A) from the bias stored in STATE (same estimator as `FES_from_State.py` of the paper repositories; checked to agree with it to all printed digits).
- `deltaF_rew`: ΔF from reweighting the COLVAR samples collected so far with `exp(bias/kBT)`.
  - Basins: `x > 0` (wolfe-quapp, mueller) and `0 < phi < 2.3` (ala2), as in the papers.
- `t_first`, `n_trans`: first passage to basin B and the number of transitions, between basin cores with hysteresis (see `opesval/systems.py`).
- ala4: number of the 8 phi-states visited over time, as in `get_ala4_basin.sh` of the Explore paper.
- Final values of `rct`, `zed`, `neff` and the number of compressed kernels.

`paper/` contains the corresponding averages from the paper repositories, for reference.
