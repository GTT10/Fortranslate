# Completion roadmap

## Definition and frozen reference

A selected PeleC responsibility is complete only when it has an independent
Fortran implementation and automated numerical evidence, or is explicitly
excluded with a reason and a usable alternative. A scaffold, test count or
clean compilation is not sufficient.

The reference is PeleC `development`
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`; recursive dependencies are pinned in
[the manifest](../references/pelec_baseline.json). Changing that baseline must
identify affected numerical signatures and acceptance criteria.

## Consolidation versus solver completion

The source-provenance repair, documentation ownership, CMake module split,
archived PR11 supersession review, and hosted serial/MPI/CVODE/clean-install
workflows are implemented. [PR #102](https://github.com/GTT10/Fortranslate/pull/102)
records the actual integration status and completed run identities; see
[current status](current_status.md) and [the layout evidence](validation/build-layout-20260909.md).
A targeted NASA7 failure-status fix and its red/green tests are recorded
[separately](validation/nasa7-conversion-status-20260909.md).

Do not redo this housekeeping as a substitute for the remaining physical
capabilities. Conversely, do not label the following capabilities implemented
merely because the repository has been consolidated.

## Ordered implementation work

| Order | Work unit | Measurable exit gate |
| --- | --- | --- |
| 1 | Validate existing physics before expanding it | Pin external cases, states, observables and tolerances; compare H2/O2 trajectories/rates with Cantera and flow fields with frozen PeleC; characterize convergence and conservation |
| 2 | Production chemistry path | Specify the target fuel mechanism and unsupported forms; validate ingestion/0D ignition; qualify stiff CFD integration with rollback and conservation before detailed-fuel claims |
| 3 | General 3D flow and mesh | Physical boundaries, then dynamic/multipatch/multilevel AMR and general EB in separate increments; require conservative regrid/reflux, changed-rank restart and convergence |
| 4 | LES and Lagrangian spray | Validate closures and particle submodels, then two-way mass/momentum/energy coupling; qualify evaporation and breakup before penetration/ignition comparisons |
| 5 | Practical release qualification | Reproducible clean install, documented supported inputs, distributed I/O/scaling measurements, physical benchmarks and owner-resolved licensing |

These are dependencies, not permission to bundle several physical models into
one PR. Existing regular 3D and static two-level AMR foundations must not be
reimplemented. Their restrictions remain explicit in current status.

## Contract for every increment

Before modifying code, specify the failure/capability, supported dimensions,
physics/parallel mode and a measurable exit test. Include a focused unit test,
an application/conservation or external-parity test, and a precise scope
statement. Separate numerical changes from structural refactoring.

For comparisons, record reference revision, mechanism/source hash, initial
state, observable, tolerance, comparison time and mesh/time-step convergence.
Do not loosen tolerances or regenerate golden outputs merely to obtain green
CI. Failed, interrupted, missing and completed checks are distinct states.
New work starts from the integrated main; a deliberate PR dependency must
state its base and merge order rather than recreating an accidental stack.

## Historical detail

The old milestone-by-milestone roadmap remains at the
[immutable checkpoint](https://github.com/GTT10/Fortranslate/blob/509d07929d52d9e570a6af1f6ec7e5836bd4fe5d/docs/completion_roadmap.md).
Detailed rationale remains in [design decisions](design_decisions/), the
[historical implementation inventory](implementation_status.md) and dated
[validation records](validation/). [Issue #101](https://github.com/GTT10/Fortranslate/issues/101)
tracks completion and actual tooling limitations. The unexecuted ancestor-only
branch-deletion audit is an administrative remainder, not a numerical feature.
