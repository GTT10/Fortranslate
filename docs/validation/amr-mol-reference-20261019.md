# Regular AMR MOL reference comparison, 2026-10-19

Recursive, one-iteration MOL agrees with the unchanged pinned PeleC HIT
application through three 100 ns steps on a spatially varying two-level grid.
Recursive, two-iteration nonreacting MOL fails the same criteria. This advances
C2 for one frozen hydro operator case; it does not close the AMR or reacting
trajectory gate. No numerical kernel or acceptance tolerance was changed.

## Implementation and scope

The former comparison script printed differences without an evolution gate.
It now runs or reuses SHA-bound native reference plotfiles and rejects wrong
physical times, hierarchy, spacing, nonfinite fields, field errors and composite
mass discrepancies. The production-integrator probe now exports physical time
and has an explicit CMake target. Inputs, pinned source hashes and acceptance
criteria are under `cases/amr_mol_reference` and
`references/pelec_amr_mol.json`.

The periodic domain is 2π cm per side, with 64 coarse cells and 256 fine cells
covering half the domain. Comparisons use 288 exposed cells. All initial cell
states are imported without projection or clipping. Chemistry, transport and
EB are absent. This is a serial operator comparison, not initializer, MPI,
regridding, convergence or physical validation.

## Results

| Recursive case | Largest scaled field error over 0–300 ns | Reference mass change at 300 ns | PeleF mass change | Verdict |
| --- | --- | --- | --- | --- |
| One MOL iteration | `1.136e-15` | `-1.858e-16` | `0` at recorded precision | Pass |
| Two MOL iterations, chemistry off | `3.397e-4` | `-1.055e-6` | `0` at recorded precision | Fail |

Initial fields agree exactly. The trajectory changes by about 1.30% of its
initial field scales, so this is not a uniform-state or zero-evolution pass.
The two-iteration difference already reaches `1.135e-4` in scaled x momentum
at 100 ns. All four time slices and all failed criteria are preserved.

Pinned `Source/Advance.cpp:99–100,133–134,170–180` registers half the old MOL
flux, suppresses the new contribution when `mol_iters > 1`, and performs the
later half registration only inside `do_react`. With chemistry off, the
accepted trapezoidal state therefore has a different coarse/fine correction.
PeleF's existing `src/amr/amr_patch_tree_split_3d_mod.F90:401–410` refluxes the
accepted state with both physical halves and reconciles the old-half register
on a private copy. The measured periodic mass difference is consistent with
this explicit policy difference. Fresh reference fields at 0, 100 and 300 ns
also exactly reproduce the retained older reference fields; reference drift
does not explain the discrepancy.

Matching this reference branch by removing the physical correction would
sacrifice conservation. The failed case remains in scope and its gate remains
red. Activation/convergence evidence for the conservation-preserving deviation
and reacting spatial-gradient AMR comparison remain open.

## Identity, commands and checks

PeleF branch `main`, base `0a0e00ebe9138145fbb0e0e745665d0d768df009`, with
pre-existing uncommitted work. `build/amr-mol-reference-20261019/` contains
`source-identity.json`, source copies, `working-tree.patch`, `git-status.txt`,
both build caches and logs. The pinned PeleC revision is
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`; its 123 relevant source hashes were
checked against the official archive. GNU Fortran/C++ 13.3.0 and CMake 3.28.3
were used: PeleF Debug, PeleC Release. Both applications ran directly in serial;
this provides no MPI or performance evidence. The surrounding PeleF build has
Open MPI and SUNDIALS enabled, but the hydro probe does not link either runtime.

Final checker records, relative to the artifact directory:

- Passed: `runs/matched_amr_mol_7kw4hoqg/checks.json`.
- Failed: `runs/matched_amr_mol_d7c74be7/checks.json`.
- Fresh native references: `runs/matched_amr_mol_q0flzfz4/reference/` (one
  iteration), `runs/matched_amr_mol_lbtv01sm/reference/` (two iterations).

The final checks reused the captured references after extending failure reports
to include every time slice. In total, two fresh native reference simulations
and four Fortran probe simulations completed: two distinct cases, one passing
and one failing. An initial invocation with a wrong probe path failed before
starting a simulation. Each `checks.json` contains exact executed commands,
binary/checker hashes, elapsed times and return codes. No common-clock cases,
full CTest matrix or installed-binary checks were run in this increment.

Completed static checks (logs beside the runs):

```sh
build/cantera-venv/bin/python tools/check_mechanism_provenance.py --bundle mechanisms/lidryer.json --source mechanisms/lidryer.yaml --regenerate
python3 tools/check_project_contract.py --source .
python3 -m py_compile tools/check_amr_tree_mol_reference.py
git diff --check
```

All four passed. The provenance check covers the selected build's pinned
mechanism; the hydro comparison itself uses a gamma-law species.


## 2026-10-28 current-source four-case audit

The production probe was rebuilt from the current dirty `main` at base
`0a0e00ebe9138145fbb0e0e745665d0d768df009`, using the existing GNU Fortran
13.3.0 configuration. All four unchanged native cases were freshly run after
checking the pinned source contract. Both solvers completed every case through
300 ns. No numerical source, input gate or reference was changed.

| Clock / MOL iterations | Maximum scaled field difference | Native final relative mass change | PeleF final relative mass change | Original gates |
| --- | --- | --- | --- | --- |
| Common / 1 | 7.5733e-16 | -1.8575e-16 | 0 at CSV precision | Pass |
| Common / 2 | 3.39190e-4 | -1.04287e-6 | 0 at CSV precision | Fail, 24 violations |
| Recursive / 1 | 1.13557e-15 | -1.8575e-16 | 0 at CSV precision | Pass |
| Recursive / 2 | 3.39739e-4 | -1.05462e-6 | 0 at CSV precision | Fail, 24 violations |

Every initial field matches exactly. The two-iteration discrepancy therefore
also occurs without recursive clocks; it cannot be attributed solely to
subcycled temporal interpolation. The fresh evidence remains consistent with
the old-half registration policy described above. Source inspection confirms
that pinned `Advance.cpp` suppresses new-time registration for `mol_iters>1`
and completes it only in `do_react`, whereas PeleF refluxes the accepted
nonreacting state with both physical halves. This is an explicit conflict
between this frozen branch's trajectory and the conservation requirement;
neither requirement is silently dropped and the case remains failed.

Artifacts: `build/amr-mol-current-matrix-20261028/`. `driver.py` preserves each
case failure and continues to the remaining cases. Its process exit zero means
the matrix finished, not that all cases passed: `all_gates_passed=false`.
`record.json` records exact commands, source hashes and timings;
`build-record.json` records the build command, compiler/cache/binary identities.
`verify.py` checks source stability, every command return code, native terminal
markers, all four field times, initial equality and the captured reference SHA
manifests. `verified.json` records the two passes and two failures separately.
The notification worker was3366033; completion was inspected directly, so no
long-job wait remains. This is serial gamma-law operator evidence, not an MPI,
EB, installed-release, convergence or physical-validation result.

Provenance and project-contract checks pass. Python syntax and diff whitespace
checks also pass. Further work must retain both the native mismatch and the
conservation invariant explicitly; the successful MOL1 branch is not a smaller
replacement for the required MOL2 comparison.


## 2026-10-28 isolated native registration control

An isolated copy of native `Advance.cpp` changes exactly one expression:
`reflux_factor = mol_iters > 1 ? 0 : 0.5;` becomes
`reflux_factor = (mol_iters > 1 && do_react) ? 0 : 0.5;`.
Only this object is recompiled and substituted into a separate executable.
The pinned source and executable remain hash-identical. This diagnostic
therefore restores the new-time half registration when chemistry is off;
it changes neither the state integration formula nor a PeleF numerical source.

All four diagnostic cases finish three steps at300 ns. With common clocks,
the maximum scaled difference from PeleF is7.5734e-16 for both iteration
settings; with recursive clocks it is1.1356e-15. Both MOL2 maxima previously
were about3.4e-4. The control's final mass drift is-1.8575e-16 with common
clocks and0 at exported precision with recursive clocks. MOL1 fields remain
exactly equal to the original fixed reference at all four exported times;
MOL2 initial fields remain equal. Scales in these diagnostic comparisons are
still the original fixed-reference field maxima, not variant-defined gates.

This single-expression intervention removes both the observed trajectory
difference and mass drift in this frozen gamma-law case, identifying the
registration condition as sufficient to explain the discrepancy. It does not
qualify the variant as a replacement reference, establish general AMR behavior,
or close the original MOL2 failures. Matching that original branch and retaining
the current periodic conservation requirement remain incompatible for this case.
No silent exclusion, tolerance relaxation or conservation-breaking compatibility
mode has been added.

Artifacts are `build/amr-mol-registration-control-20261028/`: `driver.py`,
`record.json`, the isolated source/object/executable, four case directories,
58 build/link input hashes in `link-inputs-retry.json`, and `verify.py` with
`verified.json`. Verification checks all22 compile/link/run/export commands,
16 exported field comparisons, terminal markers, control identity, and baseline
preservation. The first attempt stopped before compilation because the driver
misspelled CMake's object path (`__` segments); its failure is retained in
`notification/`. The corrected attempt is recorded in `notification-retry/`
(worker3367320). Python syntax, LiDryer provenance, project contract and diff
whitespace checks pass. These checks are diagnostic evidence, not new release
or physical-validation passes.


## 2026-10-28 preserve complete failed comparison matrices

`tools/check_amr_tree_mol_reference.py --keep-going` now completes all requested
clock/iteration comparisons after numerical field/inventory gate violations.
Completed cases have explicit passed/failed status. `matrix_complete` records
whether all requested comparisons completed, separately from overall status;
any gate violation still produces failed status and a nonzero process exit.
Provenance errors, failed commands and malformed data remain fatal. Default
fail-fast behavior and every original acceptance threshold are retained.

Four integration checks use the rebuilt production probe and the four original
native reference captures from the current-source audit, merged only after
checking common source-contract/binary identity and every artifact hash:

- Complete matrix: four cases compared, two failed, exit1, matrix complete.
- Default behavior: stops after the first failed case (two cases), exit1.
- Single passing recursive MOL1 case: exit0, matrix complete.
- Injected failing probe executable: exit1, no completed cases, matrix incomplete.

Every per-step error and violation in the complete matrix exactly matches the
preceding separate-case audit. The control-flow improvement neither masks the
MOL2 discrepancy nor treats completing a matrix as qualification. Commands,
reports and hashes are under `build/amr-mol-keep-going-20261028/validation.json`;
`validate.py` reproduces these integration checks. The failure-probe check uses
`/usr/bin/false` and is a deliberate harness failure, not a solver failure.
Python syntax, LiDryer provenance, project contract and diff whitespace checks
pass. No production numerical code changed.
