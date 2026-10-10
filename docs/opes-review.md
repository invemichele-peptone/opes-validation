# OPES module review: bugs and speed-ups

This note reviews the OPES module of PLUMED (`src/opes`, master c3506e974, which matches
v2.9 for OPES apart from documentation and keyword API changes). It looks for bugs and for
speed-up opportunities. The neighbor list has its own note, [nlist.md](nlist.md), and the
kernel cutoff and merging are covered in [kernel-cutoff.md](kernel-cutoff.md) and
[kernel-merging.md](kernel-merging.md).

How it was done:
- OPES_METAD and OPES_METAD_EXPLORE (`OPESmetad.cpp`) were reviewed by hand.
- OPES_EXPANDED and the expansion CVs (ECVs) were reviewed by a separate pass.
- Every finding listed as **reproduced** below was run and its effect observed, with a
  PLUMED built with or without MPI (OpenMPI 4.1.6). Findings marked **inspection** were
  checked by reading the code only.

**Summary.**

| id | where | what | effect | evidence |
|---|---|---|---|---|
| M1 | OPES_METAD(_EXPLORE), restart from KERNELS with adaptive σ | σ0 is reset to the first kernel's already shrunk σ | after the restart, new kernels get the wrong weight relative to old ones: ×0.79 for explore with 2 CVs, ×0.67 with 3 | reproduced |
| M2 | OPES_METAD and OPES_EXPANDED, `WALKERS_MPI` with STATE files | walkers other than the first never clear the STATE header fields | the cost of writing STATE grows quadratically with the number of writes: 154 s vs 2 s for 2000 writes | reproduced |
| E1 | ECV_CUSTOM, 2 or more ARGs, no `ADD_P0` | wrong index for the ARGs not in the ensemble | every ensemble but the first gets β·ARG0 added: wrong ΔF, bias and forces | reproduced |
| E2 | OPES_EXPANDED, observation phase, more than one MPI process per replica | STATE written before initialisation | crash in MPI_Allgatherv | reproduced |
| E3 | ECV_CUSTOM restart, `_` in ARG labels | the code splits labels on `_` | restart from STATE or DELTAFS aborts | reproduced |
| E4 | OPES_EXPANDED, `CALC_WORK` with `WALKERS_MPI` | ECVs are left at the last walker's CVs | wrong `work` on every walker but the last; bias and forces are fine | reproduced (after the first version of this note) |
| E5 | ECV_MULTITHERMAL_MULTIBARIC, `CUT_CORNER` | `temp_high-pres_low` instead of `pres_high-pres_low` | the cut line is shifted; with non-default units it can exclude no points at all | inspection |
| E6 | ECV_UMBRELLAS_LINE/FILE with several CVs | `BARRIER` caps each CV separately | the initial bias can reach (number of CVs)×BARRIER | inspection, intent unclear |
| M3 | OPES_METAD, `getMergeableKernel` with NLIST | early exit reads the shared minimum instead of the thread's own | data race under OpenMP; result unaffected, weaker early exit | inspection |
| M4 | OPES_METAD, Z update | `ks*ks` overflows `unsigned` from 65,536 kernels | at some kernel counts Z is recomputed with N² evaluations (tens of seconds); result unaffected | reproduced (arithmetic) |
| M5 | OPES_METAD, adaptive σ | `update()` returns before the STATE write during the first ADAPTIVE_SIGMA_STRIDE steps | no STATE if a checkpoint falls in that window | inspection |
| S1 | OPES_METAD kernel sum | per-kernel early exit and data layout | a two-pass version is 2.7–3.7× faster on the kernel sum (ala4), with the same results to 4·10⁻¹⁶; in MD, ala4 runs 2.0× faster (OPES_METAD) and 1.4× faster (EXPLORE) | micro-benchmark, MD |

M2 and E1–E6 come from the review of OPES_EXPANDED and the ECVs; M2 was then found and
measured in OPES_METAD as well. Line numbers refer to master c3506e974.

**Status: fixed on the plumed2 branch `opes-force-shifted-kernel`**, each with a regtest
that fails without the fix:

| commit | what |
|---|---|
| f4c856b | M1–M4 in OPES_METAD; new regtest `rt-opes_metad-kernels-restart` |
| d439427 | E1–E5, and M2 in OPES_EXPANDED; new regtests `rt-ecv_custom-nop0` (2 MPI processes), `rt-ecv_custom-nop0-restart`, `rt-opes_expanded-mw-work` |
| c0b558a | exact recursive merges with `NLIST_REBUILD=EXACT` (see [nlist.md](nlist.md)) |
| 4c1f694 | S1, the two-pass kernel sum (section 3) |
| b880f90 | faster neighbor list rebuild, with the block code of the sum (see [nlist.md](nlist.md), section 5) |

- E4 is now also reproduced: in a two-walker run, swapping the trajectories should swap
  the walkers' `work`. Before the fix, walker 0 got 16.5 instead of 2.36 on the same
  trajectory; after it, the values swap exactly.
- The E5 fix changes the references of the four regtests that use CUT_CORNER. The new
  outputs keep exactly the grid points on or above the cut line, recomputed
  independently; the old ones kept 1–6 points below it.
- Not changed: E6, whose intent is the author's call, and M5, which would need STATE
  files without kernels.

## 1. OPES_METAD and OPES_METAD_EXPLORE

### M1. Restart from KERNELS resets σ0 with adaptive σ

`OPESmetad.cpp:727-729`: when restarting from a KERNELS file (not a STATE file) with
`SIGMA=ADAPTIVE`, σ0 is taken from the first kernel in the file:

```c++
if(counter_==(1+walker_rank_) && adaptive_sigma_) {
  sigma0_=sigma;
}
```

But the first kernel's σ was already shrunk by the rescaling factor
$s_1=(N_1(d+2)/4)^{-1/(4+d)}$, where $d$ is the number of CVs:
- $N_1=2$ (the counter) for OPES_METAD_EXPLORE;
- $N_1\approx 1$ ($n_{\rm eff}$) for OPES_METAD.

So σ0 comes back too small by $s_1$. σ0 enters only the kernel heights,
$h=w\prod_i\sigma_{0,i}/\sigma_i$. New kernels therefore carry a weight $s_1^d$ relative
to the kernels read from the file:

| weight of new kernels vs old ones | 1 CV | 2 CVs | 3 CVs |
|---|---|---|---|
| OPES_METAD_EXPLORE | 0.92 | 0.79 | 0.67 |
| OPES_METAD | 1.06 | 1.00 | 0.91 |

Reproduced on ala2 OPES_METAD_EXPLORE, quick tier replica 0. Its STATE file has
σ0 = (1.8056, 2.8593). After a restart from its KERNELS file, the STATE written at
start-up has σ0 = (1.6086, 2.5474). That is a factor 0.891 = 2^(−1/6) per CV, as
predicted. Z is restored correctly (1.12628 vs 1.12630).

- PLUMED already warns that a KERNELS restart "might be approximate". But this error does
  not shrink with the length of the first run, and it is easy to avoid.
- The ratio of the first kernel's height to its weight is exactly $s_1^{-d}$, and both
  are in the file (`height`, `logweight`). So $\sigma_0=\sigma_1 (h_1 e^{-\log w_1})^{1/d}$
  recovers σ0 exactly, unless SIGMA_MIN clipped σ1.
  - For OPES_METAD_EXPLORE the stored height is already $s_1^{-d}$.
- With multiple walkers the same applies to line `1+walker_rank_`.

### M2. STATE files with multiple walkers: quadratic cost

`OPESmetad.cpp:1586-1590` (and `OPESexpanded.cpp:855-858`):

```c++
if(storeOldStates_) {
  stateOfile_.clearFields();
} else if(walker_rank_==0) {
  stateOfile_.rewind();   // also clears the fields
}
```

Walkers other than the first write their STATE to `/dev/null` and neither rewind nor
clear the fields. Every write then calls `addConstantField` about 9 + 2×(number of
periodic CVs) more times. `OFile::addConstantField` appends without checking for
duplicates. Since `zed` changes at every write, the whole header with all accumulated
fields is printed again, so write $n$ prints O($n$) header lines.

Measured with `plumed driver --multi 2` on ala2 trajectories (OPES_METAD, PACE=1), with
and without `WALKERS_MPI`:

| STATE writes | 2 walkers sharing the bias | 2 independent replicas |
|---|---|---|
| 500 | 9.1 s | 0.8 s |
| 1000 | 34.0 s | 1.2 s |
| 2000 | 153.7 s | 2.0 s |

- With the usual STATE_WSTRIDE (about 100 writes per run) this is negligible.
- With frequent writes over a long run it dominates. The memory also grows without bound.
- Fix: `else { if(walker_rank_==0) stateOfile_.rewind(); else stateOfile_.clearFields(); }`,
  in both actions.

### M3. Data race in `getMergeableKernel` with NLIST

`OPESmetad.cpp:1448`: in the NLIST branch, the early exit of the loop over CVs compares
with the shared `min_norm2`, which other threads update in the critical section. The
non-NLIST branch correctly uses the thread's own `min_norm2_omp`.
- The race is undefined behaviour, and the early exit is weaker: until a thread finishes,
  it compares with the threshold.
- The chosen kernel is still correct. A wrongly accepted partial distance is always at
  least the shared minimum already found, so it never wins the final comparison.
- Fix: `if(norm2>=min_norm2_omp)`.

### M4. Integer overflow in the choice of the Z update

`OPESmetad.cpp:1188`:

```c++
const unsigned ks=kernels_.size();
const unsigned ds=delta_kernels_.size();
const bool few_kernels=(ks*ks<(3*ks*ds+2*ds*ds*NumParallel_+100));
```

`ks*ks` wraps around in 32 bits from 65,536 kernels on. At some kernel counts the test
then becomes true, and Z is recomputed from scratch with N² kernel evaluations instead
of the incremental update.
- With one new kernel per deposition (ds = 1), the counts are 65,536–65,537, 92,682–92,683,
  113,512–113,513, …; 35 values up to 300,000. More with ds = 2–3.
- Each occurrence costs 4.3·10⁹ kernel evaluations or more, tens of seconds on one core.
  Z is the same either way (to round-off).
- Kernel counts this large are rare, but possible with many CVs or without compression.
- Fix: do the test in `double` or 64-bit integers.

### M5. No STATE during the adaptive-σ warm-up

`OPESmetad.cpp:1020-1022`: during the first ADAPTIVE_SIGMA_STRIDE steps (default
10×PACE), `update()` returns before reaching the STATE write at the end. With the default
checkpoint-only STATE writing, a GROMACS checkpoint in that window leaves no STATE file,
and a restart from it fails. This is minor: the window is short.

### Checked and found correct

- Bias, forces and the Z incremental update, including the double counting correction
  for the changed kernels.
- Recursive merging and the index bookkeeping of the neighbor list after an erase. The
  bookkeeping relies on the list being sorted, which every code path keeps.
- The MPI and OpenMP reductions and the multiple-walker exchange of kernels and weights.
  Walkers stay bitwise identical (see the parallel tests in [nlist.md](nlist.md)).
- The neighbor list paths, apart from the v2.9 inexactness discussed in
  [nlist.md](nlist.md).

## 2. OPES_EXPANDED and the expansion CVs

### E1. ECV_CUSTOM without ADD_P0 gives wrong ensembles

`ECVcustom.cpp:164-173`: `getIndex_k()` leaves `index_k[k][j]=0` for every ARG j that does
not belong to ensemble k. Index 0 should point to a zero entry, and it does when `ADD_P0`
is used. Without `ADD_P0`, entry 0 of ARG 0 holds β·ARG0. So ensemble k ≥ 1 gets
$\Delta u_k=\beta(\mathrm{ARG}_0+\mathrm{ARG}_k)$ instead of $\beta\,\mathrm{ARG}_k$.

Reproduced with `ECV_CUSTOM ARG=t1,t2 TEMP=300`, OPES_EXPANDED with `PACE=1
OBSERVATION_STEPS=1`, where t2 = 1.75·t1 by construction:

| first DELTAFS line | DeltaF_t1 | DeltaF_t2 | expected DeltaF_t2 |
|---|---|---|---|
| without `ADD_P0` | −12.402 | **−34.106** | −21.704 |
| with `ADD_P0` | −12.402 | −21.704 | −21.704 |

- The error is silent. Both regtests (`rt-ecv_custom`, `rt-ecv_custom-restart`) use
  `ADD_P0`, so neither catches it.
- Fix: drop the override; the default `index_k[k][j]=k` is correct for ECV_CUSTOM,
  because `ECVs_[j][k]` is zero for k ≠ j + P0.

### E2. STATE written during the observation phase crashes with MPI

`OPESexpanded.cpp:853-870`: during `OBSERVATION_STEPS`, the DeltaF arrays and the MPI
counts are still empty. A STATE write on a step that is not a multiple of PACE reaches
`comm.Allgatherv(deltaF_,all_deltaF_,&all_size_[0],&disp_[0])` with empty vectors.

- This happens with checkpoint-only STATE writing when a checkpoint falls in the
  observation phase, or with a STATE_WSTRIDE that is not a multiple of PACE. The check at
  line 241 only requires STATE_WSTRIDE ≥ PACE, although its message says "multiple".
- Reproduced: `PACE=2 OBSERVATION_STEPS=10 STATE_WFILE=State.data STATE_WSTRIDE=3`
  - runs with 1 MPI process, writing an empty state;
  - with `mpirun -np 2` it aborts with "An error occurred in MPI_Allgatherv".
- Fix: skip the STATE write until the action is initialised, and require STATE_WSTRIDE
  to be a multiple of PACE.

### E3. ECV_CUSTOM restart fails with `_` in ARG labels

The DeltaF column names are built from the ARG labels (`DeltaF_t_1_NaN`), and restarts
parse them by splitting on `_` (`OPESexpanded.cpp:337-358`, `871-878`).
- Reproduced: the same input restarts fine with labels `t1,t2`.
- With labels `t_1,t_2`, restarts from both STATE and DELTAFS abort with "RESTART -
  mismatch in lambda values of ECV_CUSTOM".
- `_` is common in labels (`dU_1`), and the failure is loud, so this is a usability bug.

### E4. `CALC_WORK` with `WALKERS_MPI`

`OPESexpanded.cpp:627-676`: with multiple walkers, the DeltaF update loops over the
walkers and calls `calculateECVs` with each walker's CVs. When the work is computed
afterwards, the ECVs hold the last walker's values. The new bias is then evaluated at the
last walker's position, but the current bias subtracted is this walker's own.
- Only the `work` output is affected; the next `calculate()` recomputes the ECVs.
- The regtest `rt-opes_expanded-mw` gives both walkers the same trajectory, so it cannot
  see this.
- Fix: recompute this walker's ECVs before the work block.

### E5. CUT_CORNER subtracts a pressure from a temperature

`ECVmultiThermalBaric.cpp:299`: `small_value=(temp_high-pres_low)/1e4`, clearly meant to
be `(pres_high-pres_low)/1e4`. The tolerance shifts the whole cut line.
- In default units the shift is small.
- In other units it can be large; it then excludes no points, or too many. For example,
  with `UNITS LENGTH=A ENERGY=kcal/mol` the review pass estimated a shift of a few kbar
  for the documentation example. That figure was not run here.

### E6. Umbrella BARRIER applied per CV

`ECVumbrellasLine.cpp:269`, `ECVumbrellasFile.cpp:254`: the initial ECVs are capped at
BARRIER/kT for each CV j separately, while Δu_k sums over the CVs. With 2 CVs, the initial
bias can reach 2×BARRIER, as in the 2D example of the documentation. Whether this is
intended is for the author to say; the documentation suggests one barrier for the whole
system.

### Not verified here

The review pass also reported the following; they were not checked:
- an off-by-one in the automatic number of steps (`ExpansionCVs.cpp:202-206`);
- DELTAFS restarts losing up to PRINT_STRIDE − 1 updates of the counter, although the
  documentation calls them equivalent to STATE restarts;
- several minor input-validation issues (PACE=0, STORE_STATES files as STATE_RFILE);
- three MPI collectives per step in OPES_EXPANDED, which may cost more than they save for
  small numbers of ensembles.

## 3. Speed of the kernel sum

Without a neighbor list, OPES_METAD sums all compressed kernels at every MD step. In the
ala4 runs this is most of the run time: 265 s for 10 ns without a list, 141 s with one
([nlist.md](nlist.md)).

A micro-benchmark (`bench/kernel_sum.cpp`, run by `bench/kernel_sum.sh`) evaluates the
bias and its derivatives with real kernels along a real trajectory. The kernels are the
final compressed ones of ala4 quick-tier replica 0, and the trajectory is that run's
20,001 frames. It compares four versions of the loop on one core (`g++ -O3`, as PLUMED is
built here):
- **current**: the code as is, with each kernel holding its own centre and σ vectors,
  `Value::difference`, a division by σ, and an exit as soon as the partial distance
  exceeds the cutoff;
- **inv_sigma**: the same, multiplying by a stored 1/σ;
- **no_exit**: the same as current, checking the cutoff once after the loop over CVs;
- **two_pass**: centres and 1/σ stored per CV in contiguous arrays. A branch-free pass
  computes the squared distance to every kernel, and a second pass evaluates the exp and
  the derivatives for the kernels within the cutoff.

| ns per kernel (speed-up) | ala4 OPES_METAD, 15,119 kernels | ala4 EXPLORE, 2,385 kernels |
|---|---|---|
| current | 8.9 | 20.9 |
| inv_sigma | 7.9 (1.13×) | 18.1 (1.15×) |
| no_exit | 8.6 (1.03×) | 12.3 (1.71×) |
| two_pass | **3.3 (2.72×)** | **5.7 (3.70×)** |
| kernels within the cutoff | 0.7% | 13.2% |

- **The exp is not the bottleneck.** Most kernels are outside the cutoff. The cost is the
  per-kernel distance test, with its unpredictable exit branch and scattered memory.
- **two_pass is 2.7–3.7× faster.** The sums agree with the current code to 4·10⁻¹⁶ in
  relative terms. The only arithmetic change is the periodic wrap, `nearbyint` instead of
  `Tools::pbc`.
  - Making the distance loop vectorise as well (by rounding with a magic constant) gave
    no further gain.
- The micro-benchmark had predicted about 1.7–1.8× for ala4 without a list. That estimate
  assumed the kernel sum took 60–70% of those runs.
- The cheap changes (1/σ, no early exit) give 1.1–1.7× in the benchmark.

### Implementation and measured speed-up

Branch commit 4c1f694 implements the two-pass sum in OPES_METAD:
- **Storage:** the kernels that are summed (all of them, or those in the neighbor list)
  are also kept in arrays with one CV per contiguous array, together with 1/σ.
  - These arrays are rebuilt whenever the kernels or the list change, i.e. at most once
    per deposition or list rebuild. That costs O(N·d), against O(N) per MD step for the
    sum itself.
- **Periodic wrap:** a branch-free rounding, `nearbyint` where it is a single
  instruction (aarch64, and x86 with SSE4.1), otherwise a conversion to `int`. On this
  machine `Tools::pbc` cancels most of the gain: the two-pass sum is then only 1.1–2.0×
  faster.
  - The x86 fallback was not measured. In the benchmark on aarch64 it gives 1.9–2.9×.
- **Order of the sum:** unchanged, so a serial run sums the same terms in the same
  order.
- **Parallel split:** MPI processes take contiguous shares of the kernels, and OpenMP
  threads blocks of 64. Threads are used above 2 × threads × 64 kernels.

Checks:
- All 45 OPES regtests pass, with 1 and with 4 OpenMP threads.
- **Replays:** in the replays of [nlist.md](nlist.md), the new build agrees with the
  previous one to round-off. The largest differences are 1.2·10⁻¹² kBT in the bias and
  8·10⁻¹² in the forces. This holds serially, with 4 threads, with 4 MPI processes and
  with walkers.
- **Exact list:** with `NLIST_REBUILD=EXACT NLIST_ZED=FULL`, the replays stay bitwise
  identical to no list in serial runs, and at round-off in parallel ones.
- **MD:** all 12 runs of 10 ns with the exact list have COLVAR files identical to those
  without a list.

**MD timing** (`build/time_builds.sh`): 10 ns, one job per core, pinned to the 10
Cortex-X925 cores. The jobs of the old and new builds are shuffled together. Times are
mean ± s.d. of 3 replicas.

| wall time | old (160bdfc) | new (4c1f694) | speed-up |
|---|---|---|---|
| ala2 OPES_METAD, no list | 68.0 ± 1.9 s | 57.5 ± 0.5 s | 1.18× |
| ala2 OPES_METAD, exact list | 59.7 ± 2.0 s | 56.2 ± 0.6 s | 1.06× |
| ala2 EXPLORE, no list | 57.0 ± 0.9 s | 54.9 ± 1.2 s | 1.04× |
| ala2 EXPLORE, exact list | 57.0 ± 0.1 s | 54.4 ± 0.5 s | 1.05× |
| ala4 OPES_METAD, no list | 211.0 ± 4.8 s | 107.0 ± 1.5 s | **1.97×** |
| ala4 OPES_METAD, exact list | 103.8 ± 2.3 s | 87.1 ± 0.6 s | 1.19× |
| ala4 EXPLORE, no list | 130.2 ± 2.7 s | 91.5 ± 1.0 s | **1.42×** |
| ala4 EXPLORE, exact list | 130.2 ± 4.4 s | 91.4 ± 2.0 s | **1.42×** |

- **ala4, OPES_METAD without a list:** the run is now 2× faster, matching what the exact
  neighbor list alone gave.
- **ala4, OPES_METAD_EXPLORE:** 1.4× faster, where the neighbor list never helped.
- **ala2:** small gains, because it has few kernels (100–700) and the MD itself dominates.
- **With the exact list (ala4 OPES_METAD):** the new build adds 1.19×. The list now gives
  1.23× over the no-list run, down from 2×.
  - The "exact" variant of the new build also includes the exhaustive recursive merge
    search of c0b558a, one O(N) distance loop per deposition.
- **Why MD beats the benchmark estimate:** if the sum got 2.7× faster, the measured 2.0×
  means it took about 78% of the old ala4 OPES_METAD run. The estimate of 60–70% used
  per-kernel costs measured on an A725 core and run times measured on both core types.
- **Code placement adds up to 10% of noise between builds.** Later concurrent A/B runs
  ([nlist.md](nlist.md), section 5) showed this for the sum without a list, even with
  identical machine code: its speed depends on where the linker places the function.
  - Built with `-falign-functions=64`, 4c1f694 spends 6.60 µs per step in OPES for ala4
    OPES_METAD without a list, against 6.01 µs for the default build used here, about
    3 s more over the 107 s run.
  - So with a less lucky placement the 1.97× would be about 1.9×. The placement of the
    old build is not known either.

## 4. Reproducing

PLUMED builds: `pr1463` and `nlist` as in [nlist.md](nlist.md); `nlist-mpi` with
`MPI=opt/openmpi build/build_plumed.sh nlist-mpi <sha> /path/to/plumed2`.

```bash
# S1, kernel-sum micro-benchmark, and MD timing of the old and new builds
bench/kernel_sum.sh nlist pr1463
build/build_plumed.sh fast 4c1f694acbfbd48a3b27f71afad0cbe815905b75 /path/to/plumed2
CORES="5 6 7 8 9 15 16 17 18 19" build/time_builds.sh runs/kernel-sum-timing nlist fast
build/replay_nlist.sh fast pr1463 runs/nlist-replay-fast
uv run python build/check_replay.py runs/nlist-replay-fast runs/nlist-replay   # new vs old build

# M1, KERNELS restart: restart from the KERNELS file of an adaptive-sigma explore run, read sigma0
cat > plumed.dat << 'EOF'
RESTART
phi: TORSION ATOMS=5,7,9,15
psi: TORSION ATOMS=7,9,15,17
opes: OPES_METAD_EXPLORE ARG=phi,psi PACE=500 BARRIER=30 TEMP=300 FILE=KERNELS.in STATE_WFILE=STATE.conv STATE_WSTRIDE=1000000
EOF
cp runs/pr1463/quick/ala2/explore/000/KERNELS KERNELS.in
plumed driver --plumed plumed.dat --mf_xtc runs/pr1463/quick/ala2/explore/000/md.xtc
grep sigma0 STATE.conv runs/pr1463/quick/ala2/explore/000/STATE | head -4

# M2, STATE writes with walkers: md.0.xtc and md.1.xtc are two ala2 trajectories; time with and without WALKERS_MPI
echo 'phi: TORSION ATOMS=5,7,9,15
psi: TORSION ATOMS=7,9,15,17
opes: OPES_METAD ARG=phi,psi PACE=1 BARRIER=30 TEMP=300 FILE=/dev/null STATE_WFILE=STATE STATE_WSTRIDE=5 WALKERS_MPI' > plumed.dat
time mpirun -np 2 plumed driver --multi 2 --plumed plumed.dat --mf_xtc md.xtc

# E1, ECV_CUSTOM without ADD_P0 (alanine.xtc from regtest/opes/rt-ecv_custom)
echo 'phi: TORSION ATOMS=5,7,9,15
ene: CUSTOM PERIODIC=NO ARG=phi FUNC=50*x+100
t1: CUSTOM PERIODIC=NO ARG=ene FUNC=(300/500-1)*x
t2: CUSTOM PERIODIC=NO ARG=ene FUNC=(300/1000-1)*x
mc: ECV_CUSTOM ARG=t1,t2 TEMP=300
opes: OPES_EXPANDED ARG=mc.* PACE=1 OBSERVATION_STEPS=1 PRINT_STRIDE=1 FMT=%.6f' > plumed.dat
plumed driver --plumed plumed.dat --mf_xtc alanine.xtc && head -3 DELTAFS   # DeltaF_t2 should be 1.75*DeltaF_t1

# E2, STATE during the observation phase
echo 'phi: TORSION ATOMS=5,7,9,15
ene: CUSTOM PERIODIC=NO ARG=phi FUNC=50*x+100
t1: CUSTOM PERIODIC=NO ARG=ene FUNC=(300/500-1)*x
mc: ECV_CUSTOM ARG=t1 TEMP=300 ADD_P0
opes: OPES_EXPANDED ARG=mc.* PACE=2 OBSERVATION_STEPS=10 STATE_WFILE=State.data STATE_WSTRIDE=3' > plumed.dat
mpirun -np 2 plumed driver --plumed plumed.dat --mf_xtc alanine.xtc   # aborts in MPI_Allgatherv

# E3: as rt-ecv_custom-restart, with ARG labels t_1,t_2 instead of t1,t2
```
