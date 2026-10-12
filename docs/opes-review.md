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
| S2 | OPES_METAD, update of Z at each deposition | the changed kernels are evaluated against all kernels one at a time, as in the old sum | with the code of S1 this update is 6–8× faster, and OPES 2–9% faster with the exact neighbor list | MD |
| S3 | OPES_METAD, what is left | two vector allocations per step; Z at a KERNELS restart from all pairs of kernels, one at a time | OPES up to 3% faster; Z at the restart 4× faster, the restart itself 3–4%. A per-kernel skin for the exact list did not pay off | MD |

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
| 403d480 | S2, the update of Z with the block code of the sum (section 3) |
| e0f56f8 | S3, no allocations at each step, and Z at a KERNELS restart with the block code (section 3) |

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

### The update of Z (S2)

Z is proportional to the sum of all kernels evaluated at all kernel centres. A deposition
changes only a few kernels: it adds one, or replaces a kernel with its merge with the new
one (more with recursive merges or walkers). So at each deposition OPES_METAD updates Z
with the change of that sum: each changed kernel evaluated at every centre, and every
kernel evaluated at the changed kernel's centre.
- **Before:** this used `evaluateKernel`, one kernel at a time, as the sum did before S1:
  about four evaluations per kernel at a deposition with a merge.
- **Branch commit 403d480** uses the block code of the sum:
  - For a block of kernels, one branch-free pass computes both squared distances to a
    changed kernel: in units of its σ and of each kernel's own σ. Only the kernels within
    the cutoff are then evaluated.
  - The kernels are the summed ones. With `NLIST_ZED=FULL` they are the per-CV copy of
    all kernels kept for the list rebuild (b880f90), which now also stores the heights.
  - MPI processes take contiguous shares of the kernels and OpenMP threads blocks of 64,
    as for the sum.
  - That copy may now be refreshed at a deposition, before the next rebuild. The rebuild
    still uses the smallest σ of the kernels at the time of the rebuild, so the list
    behaves as before.

Checks:
- All 45 OPES regtests pass, with 1 and 4 threads.
- **Replays** of [nlist.md](nlist.md), serially and with 4 threads, 4 MPI processes,
  4 walkers, and 2 walkers × 2 processes × 2 threads:
  - All four actions agree with b880f90 to 2·10⁻¹¹ kBT in bias and forces, and to 10⁻¹³
    in Z. Only the order of the sums changed.
  - Kernel counts, list sizes and rebuild intervals are identical.
  - The exact list stays bitwise identical to no list in serial runs and with 4 walkers,
    and agrees with it to round-off with MPI.
- **MD** runs drift apart from those of b880f90, as any round-off change makes them do,
  e.g. 1568 vs 1578 kernels at the end of a Müller–Brown 2D run.

**MD timing**, concurrent A/B of b880f90 and 403d480 (`build/ab_nlist_cases.sh`, cases of
[nlist.md](nlist.md), section 5):
- Both builds use `-falign-functions=64`, against code-placement noise.
- They run at the same time on the two halves of the X925 cores, then swap halves.
- 6 runs per case and build, 4 for the long run.
- "update" is PLUMED's update step, whose fixed cost is about 0.067 µs; "OPES" is the
  whole OPES time, as in [nlist.md](nlist.md).

| µs per MD step, b880f90 → 403d480 | kernels | update, no list | update, exact list | OPES, no list | OPES, exact list |
|---|---|---|---|---|---|
| ala2 φ, OPES_METAD | 43 | 0.067 → 0.066 | 0.068 → 0.067 | 0.62 → 0.61 | 0.63 → 0.63 |
| ala2 φ,ψ, OPES_METAD | 768 | 0.124 → 0.106 | 0.120 → 0.104 | 1.35 → 1.33 | 0.91 → 0.89 |
| Müller–Brown 2D, OPES_METAD | 1570 | 0.131 → 0.100 | 0.126 → 0.099 | 1.89 → 1.85 | 1.04 → 1.02 |
| ala2 φ,ψ, jpcl | 3262 | 0.218 → 0.165 | 0.207 → 0.162 | 3.08 → 3.03 | 1.17 → 1.12 |
| ala2 φ,ψ, jpcl, long | 10⁴ | 0.614 → 0.437 | 0.576 → 0.416 | 8.88 → 8.68 | 2.22 → 2.01 |
| ala4 3 CVs, EXPLORE | 1563 | 0.186 → 0.122 | 0.188 → 0.124 | 3.10 → 3.02 | 3.05 → 2.98 |
| ala4 3 CVs, OPES_METAD | 8010 | 0.263 → 0.163 | 0.237 → 0.146 | 6.64 → 6.51 | 1.71 → 1.61 |

- **The update step** loses 30–55% of its cost above the fixed one.
- **OPES** gets 2–9% faster with the exact list and about 2% faster without, since both
  lose about the same absolute time. The wall time of the long jpcl run drops by 1%: 1032 →
  1022 s without a list, 695 → 686 s with it.
- **The speed-up of the exact list** grows a little, e.g. from 4.0× to 4.3× in the long
  jpcl run. These builds place the sum without a list less luckily than the 4c1f694 build
  used in [nlist.md](nlist.md), which gave 3.5× there.

**Where the rest goes.** A build with timers around the parts of a deposition (not
committed) ran the short jpcl and ala4 cases once. Means over the run, in µs per
deposition:

| | jpcl, no list | jpcl, exact list | ala4, no list | ala4, exact list |
|---|---|---|---|---|
| adding the kernel: merge search and KERNELS file | 14.8 | 8.7 | 19.7 | 9.1 |
| copy of the kernels into the per-CV arrays | 3.1 | 5.7 | 8.7 | 11.5 |
| update of Z | 5.4 | 5.4 | 8.0 | 8.1 |
| list rebuilds between two depositions (calculate step) | – | 45 | – | 34 |

- **The update of Z itself is 6–8× faster.** Its old cost is the new one, plus the time
  saved, plus the copy, which the sum used to make at the next step.
- **The copy** is not new work: before, it was made at the next step (no list) or at the
  next rebuild (exact list).
- **The merge search is now the largest part of a deposition.** It still compares one
  kernel at a time with `Value::difference`. It goes through all kernels twice per merge
  without a list, and once with the exact list, whose recursive merges search all kernels
  (c0b558a).
  - The block code could do this too, but the per-CV copies would have to be patched
    after each merge.
  - At ≤20 µs per deposition, i.e. ≤0.04 µs per step with PACE=500, this is not worth it
    for now.
- **In the long jpcl run, most of the remaining update time is writing STATE.** That input
  stores the state every 100 depositions (`STORE_STATES`). Each write takes up to 27 ms
  with 10⁴ kernels, 0.2–0.3 µs per step in total, for a 1 GB file.
- **With the exact list, the rebuilds cost more than the whole update.** They happen in
  the calculate step: every ~37 steps for jpcl, about 3 µs each.

### What else could be sped up (S3)

**Most of the remaining fixed cost is not OPES.** With a RESTRAINT of zero strength in place
of a 1D OPES_METAD (ala2, 403d480):
- The RESTRAINT costs 0.49 µs per step: 0.36 in calculate, 0.08 in apply, 0.04 in update.
  That is the cost of running any bias action in PLUMED.
- OPES_METAD costs 0.62 µs, so only 0.13–0.19 µs is specific to OPES.
- The whole PLUMED step costs about 4.5 µs in both cases, with two TORSIONs and a PRINT.
  The MD step costs 11–14 µs for ala2 in vacuum and 170–350 µs for a solvated GPU run.

**The candidates** (µs per MD step):

| what | cost now | idea | expected gain | drawbacks | status |
|---|---|---|---|---|---|
| rebuilds of the exact list, 2 CVs | ~0.6 at 10⁴ kernels (a rebuild every 18 steps), ~0.1 at 3000 | (a) a skin per rebuild, from the distance and σ of the kernels outside the list | up to half | a subtler exactness argument, more work per rebuild | tried, not adopted (below) |
| | | (b) a grid over CV space, for up to 3 CVs | most of it | much more code: periodic CVs, kernels whose σ differs up to 5×, merges that move centres | not done |
| exp of the listed kernels within the cutoff | 0.1–0.3 | a vectorised exp (glibc libmvec, SLEEF) | about half | compiler flags and libraries, changes in the last bit, portability | not done |
| STATE writes | 0.2–0.3, only when written often (the jpcl input: every 100 depositions) | format each kernel line at once, instead of field by field (~2.5 µs per line now) | several × on the writes | bypasses the field bookkeeping of PLUMED's files; the format must not change | not done |
| merge search and per-CV copies at a deposition | ≤0.04 and 0.02–0.05 | block code for the merge search, kernels stored per CV | most of it | a refactor of the kernel storage, which restart, STATE, walkers and merges use | not done |
| OPES's own per-step cost | 0.13–0.19 | no allocations, branch-free rebuild check | a few hundredths | none | done, e0f56f8 |
| Z at a KERNELS restart | 0.8 s once, at 10⁴ kernels | block code of the sum | 4× | none | done, e0f56f8 |

**(a) A skin per rebuild: not adopted.** The EXACT rule allows a displacement of
$(\sqrt{P_0}-1)\,r_c$ in units of the smallest σ of all kernels. A tighter bound uses each
kernel outside the list:
- A kernel at distance $d$ in units of its own σ gets within the cutoff only after a
  displacement of $d-r_c$ in those units. That is at least $(d-r_c)/\rho$ in units of the
  smallest σ, where $\rho$ is the largest ratio of the smallest σ to its σ over the CVs.
- The skin is the minimum of this over the kernels outside the list. It is never smaller
  than the current one, and the list stays exact.
- Without kernels nearby the list would never be rebuilt, so the skin is capped at a factor
  $f$ times the current one.

A temporary build (`build/patches/nlist-skin-scan.patch`) read $f$ from the environment.
$f$ = 1 is the current rule. Exact list, 2 runs per case (rebuild interval in steps / listed
kernels / OPES µs per step):

| case | f = 1 | f = 2 | f = 4 |
|---|---|---|---|
| ala2 φ,ψ, jpcl | 37 / 102 / 1.09 | 43 / 102 / 1.11 | 43 / 102 / 1.17 |
| Müller–Brown 2D, OPES_METAD | 49 / 137 / 0.99 | 57 / 137 / 1.01 | 87 / 138 / 1.05 |
| ala4 3 CVs, OPES_METAD | 102 / 327 / 1.57 | 117 / 327 / 1.62 | 2553 / 2214 / 4.12 |
| ala2 φ,ψ, OPES_METAD | 70 / 79 / 0.85 | 2549 / 297 / 1.08 | 5002 / 515 / 1.30 |

- **Where many kernels lie just outside the list**, they set the skin, and rebuilds are only
  ~15% rarer. Probably kernels with about the smallest σ are always among them, so the
  bound hardly improves. The extra work in each rebuild cancels the gain.
- **Where the cap is reached instead**, the list is rarely rebuilt and keeps every new
  kernel, and OPES gets slower: up to 2.6× for ala4 at $f$ = 4.
- So the rule of b880f90 stays.

**Per-step cost (e0f56f8).** `calculate()` now keeps the CV values and the derivatives in
member buffers, instead of allocating two vectors at every step. The EXACT check also wraps
periodic CVs without branches, as the sum does. Concurrent A/B against 403d480, aligned
builds, 6 runs each, OPES µs per step:

| case | no list | exact list |
|---|---|---|
| ala2 φ, OPES_METAD | 0.609 → 0.601 | 0.621 → 0.612 |
| Müller–Brown 1D, OPES_METAD | 0.644 → 0.631 | 0.634 → 0.617 |
| ala2 φ,ψ, OPES_METAD | 1.328 → 1.308 | 0.886 → 0.866 |
| ala2 φ,ψ, jpcl | 3.042 → 3.015 | 1.111 → 1.087 |
| ala4 3 CVs, OPES_METAD | 6.515 → 6.509 | 1.597 → 1.581 |

- OPES is up to 0.03 µs per step faster, up to 3% of its time. Wall times do not change
  measurably.
- Most of what OPES itself still costs at each step is probably the sum over the listed
  kernels. In 1D about 6 of the ~10 listed kernels are within the cutoff, and each needs
  an exp.

**Z at a KERNELS restart (e0f56f8).** Without a STATE file, a restart reads all deposited
kernels from the KERNELS file, merges them again and computes Z from scratch. Z is now
computed with the kernel sum used for the bias, instead of all pairs of kernels one at a
time.
- Test: restart from the KERNELS file of a long jpcl run, with 10⁵ deposited kernels
  compressed to 10⁴, on one Cortex-A725 core. Each build also ran with `NO_ZED`, which
  skips Z, to time Z alone.
- Z takes 0.2 s instead of 0.8 s, and agrees to all 15 printed digits.
- The restart as a whole takes 13 s and is only 3–4% faster. Most of it is reading the
  file and merging every kernel again. A STATE file avoids both.

Checks of e0f56f8:
- All 45 OPES regtests pass, with 1 and 4 threads.
- The replays are bitwise identical to 403d480 serially, with 4 MPI processes and with
  4 walkers.
- With OpenMP they differ by up to 9·10⁻¹². The order of the OpenMP reductions changes
  from run to run.
- Kernel counts, list sizes and rebuild intervals are identical.

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

# S2, update of Z: replays as for S1 (also with THREADS, NP, WALKERS), and the concurrent A/B
build/build_plumed.sh zfast 403d480b42b21021bc00389119a7abbef2501265 /path/to/plumed2
build/replay_nlist.sh zfast pr1463 runs/nlist-replay-zfast
uv run python build/check_replay.py runs/nlist-replay-zfast runs/nlist-replay-fastnl   # vs b880f90
CXXFLAGS="-O3 -falign-functions=64" build/build_plumed.sh fastnl-a64 b880f90cf2d59fd1268c8f8c2ad90ea7d09a5013 /path/to/plumed2
CXXFLAGS="-O3 -falign-functions=64" build/build_plumed.sh zfast-a64 403d480b42b21021bc00389119a7abbef2501265 /path/to/plumed2
build/ab_nlist_cases.sh runs/zed-ab fastnl-a64 zfast-a64 ala2-1d-opes mueller2d-opes ala2-2d-opes ala2-2d-jpcl ala4-3d-opes ala4-3d-explore
SCALE=10 REPS=2 build/ab_nlist_cases.sh runs/zed-ab-long fastnl-a64 zfast-a64 ala2-2d-jpcl
uv run python build/summarize_nlist_cases.py --pool runs/zed-ab/zfast-a64-1 runs/zed-ab/zfast-a64-2

# S3: per-step cost and KERNELS restart (e0f56f8), and the scan of the per-kernel skin
build/build_plumed.sh zfast2 e0f56f8ca55d8a8df84fcbc9b84142f53a00ef5d /path/to/plumed2
CXXFLAGS="-O3 -falign-functions=64" build/build_plumed.sh zfast2-a64 e0f56f8ca55d8a8df84fcbc9b84142f53a00ef5d /path/to/plumed2
build/replay_nlist.sh zfast2 pr1463 runs/nlist-replay-zfast2
uv run python build/check_replay.py runs/nlist-replay-zfast2 runs/nlist-replay-zfast
build/ab_nlist_cases.sh runs/floor-ab zfast-a64 zfast2-a64 ala2-1d-opes mueller1d-opes ala2-2d-opes mueller2d-opes ala2-2d-jpcl ala4-3d-opes
# the skin scan: commit the patch on a branch from e0f56f8, build it as skinscan-a64, then for f in 1 2 4 8
git -C /path/to/plumed2 apply /path/to/opes-validation/build/patches/nlist-skin-scan.patch
OPES_SKIN_FACTOR=$f VARIANTS=exact REPS=2 CORES="5 6 7 8 9" build/time_nlist_cases.sh runs/skinscan/f$f skinscan-a64 \
  ala2-1d-opes ala2-2d-opes mueller2d-opes ala2-2d-jpcl mueller2d-explore ala4-3d-opes ala4-3d-explore
# KERNELS restart: KERNELS.in from runs/zed-ab-long/zfast-a64-1/ala2-2d-jpcl-exact-0, timed with and without NO_ZED
echo 'RESTART
phi: TORSION ATOMS=5,7,9,15
psi: TORSION ATOMS=7,9,15,17
opes: OPES_METAD ARG=phi,psi TEMP=300 PACE=500 SIGMA=0.15,0.15 BARRIER=50 BIASFACTOR=10 CUTOFF_SCHEME=PR NLIST NLIST_REBUILD=EXACT NLIST_ZED=FULL FILE=KERNELS.in STATE_WFILE=STATE STATE_WSTRIDE=-1
PRINT FMT=%.15g STRIDE=1 FILE=COLVAR ARG=phi,psi,opes.*' > plumed.dat
taskset -c 0 opt/plumed/zfast2-a64/bin/plumed driver --plumed plumed.dat --mf_xtc runs/pr1463/quick/ala2/opes/000/md.xtc

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
