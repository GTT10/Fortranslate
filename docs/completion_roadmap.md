# Completion roadmap

## Definition and frozen reference

A selected PeleC responsibility is complete only when it has an independent
Fortran implementation and automated numerical evidence, or is explicitly
excluded with a reason and a usable alternative. A scaffold, count of tests,
or clean compilation is not sufficient.

The reference is PeleC `development`
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`; recursive dependencies are pinned in
[the manifest](../references/pelec_baseline.json). Baseline changes must state
which numerical signatures and acceptance criteria change.

## Ordered work

| Priority | Implementation unit | Exit gate |
| --- | --- | --- |
| 0 | Complete PR #102 checkpoint integration: source provenance, CI/install contracts, evidence reconciliation | Repaired source identified; required serial/MPI/selected/CVODE configurations and separate install gates completed or explicitly left open; no failing required CI; historical results not relabeled |
| 1 | Finish development-system consolidation | Keep README/current status/roadmap distinct; split large CMake target/test registrations in behavior-preserving increments; demonstrate identical target/test manifests and numerical outputs before and after |
| 2 | Validate existing physics before expanding it | Pin external cases, initial states, observables and tolerances; compare H2/O2 trajectories/rates with Cantera and relevant flow fields with frozen PeleC; characterize convergence, conservation and limits |
| 3 | Production chemistry path | Select a target fuel mechanism; enumerate unsupported reaction/thermo/transport forms; validate ingestion and 0D ignition first; qualify a stiff integrator in CFD with rollback and conservation before detailed-fuel claims |
| 4 | Generalize 3D flow and mesh capabilities | Implement physical boundaries, then dynamic/multipatch/multilevel AMR and general EB in isolated increments; require conservative regrid/reflux, rank changes, restart and convergence tests |
| 5 | LES and Lagrangian spray | Validate each closure and particle submodel independently, then two-way mass/momentum/energy coupling; test evaporation and breakup before spray penetration/ignition comparisons |
| 6 | Practical release qualification | Reproducible clean install, documented supported inputs, scalability/I/O measurements, physical benchmarks, and owner-resolved licensing |

Priorities are dependencies, not permission to bundle several workstreams into
one PR. Do not reimplement the already existing regular 3D or static two-level
AMR foundations. See [current status](current_status.md) for their exact limits.

## Contract for every new increment

Before code changes, specify the failure or capability, supported dimensions,
physics and parallel mode, and one measurable exit test. Include a focused
unit test, an application/conservation or external-parity test, and an exact
scope statement. Separate numerical changes from structural refactoring.

For comparisons, record the reference revision, mechanism/source hash, initial
state, observable, tolerance, comparison time and mesh/time-step convergence.
Do not repair failures by loosening tolerances or regenerating golden results
without explaining the numerical reason. Mark interrupted and unavailable
checks as incomplete, not passed.

## Historical detail

The previous milestone-by-milestone roadmap is preserved at the
[immutable checkpoint](https://github.com/GTT10/Fortranslate/blob/509d07929d52d9e570a6af1f6ec7e5836bd4fe5d/docs/completion_roadmap.md).
Detailed rationale remains in [design decisions](design_decisions/), the
[implementation inventory](implementation_status.md), and dated
[validation records](validation/). New task ownership and completion tracking
belong to [issue #101](https://github.com/GTT10/Fortranslate/issues/101), not
repeated copies of the same milestone narrative in every document.
