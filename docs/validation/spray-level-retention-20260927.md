# Retained physical level stages — 2026-09-27

Implementations `2c37f882` (routing) and `21c85621` (stage), based on `b8db349a`,
local `main`, version 0.246.0.
Pinned PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Implementation

`redistribute_mpi_amr_level_particles` keeps live physical parcels on the
specified level while permitting motion within its grown patch boxes. It
shares validated exchange with ghost routing, but collectively rejects a live
parcel outside every grown box instead of deleting it. Root-level routing
uses capped root lookup, including covered cells. Dead particles are removed
before positional validation. Parcel, component and node arrays publish together.

`advance_mpi_amr_level_halo` connects this routing after the existing frozen-gas
kick, drift and same-level source transfer. A routing failure leaves particles,
components, nodes, valid-cell source, halo remainder and boundary/wall ledgers
unchanged. The physical and existing ghost stage APIs share one private stage
implementation; ghost deletion semantics remain unchanged. Gas halos stay
read-only. Source density is the time-integrated reaction for the supplied kick.
No coarse source feedback is enabled by this stage.

The caller must keep motion inside its level grow region, or handle rejection
by reducing the interval. Final synchronized collection changes level ownership.
This prepared-halo stage does not schedule breakup, child creation, films,
external acceleration or SDC population lifecycles. Those responsibilities and
recursive spray application integration remain open.

## Completed checks

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6. Four focused tests passed, with no failures/skips:

- `spray_mpi_amr_tree_2`, `_4`: 0.48/0.50 s; CTest total 1.03 s.
- `spray_mpi_amr_step_2`, `_4`: 0.82/0.78 s; CTest total 1.65 s.

The routing fixture verifies every physical ID remains exactly once after
crossing a fine valid edge into its grow region. Moving one live parcel beyond
all grown fine boxes rejects collectively and preserves payload/coordinates.
Root routing retains all live particles over refinement and removes a dead
particle before checking its out-of-domain position.

The stage fixture exercises drag, heating, evaporation and full drift with the
kicked velocity. A final exchange capacity failure preserves particle motion,
components, node assignments, valid sources, remainders and ledgers. Existing
ghost-stage success checks still execute through the shared implementation.
An actual grow-region escape is checked at routing level, not independently
through the full kick/drift stage.

Optional EB physical-stage inputs compile but were not independently run here.
No complete recursive particle lifecycle, full matrix, installed audit, external
solver numerical parity or physical validation was performed. Provenance,
Python syntax, project contract and whitespace checks passed.

## Reproduction

Artifacts: `build/spray-level-retention-checks/`, including build and stage-build
logs, provenance, contract, CTest summaries and copied detailed logs.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_tree --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_tree_[24]$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_step --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_step_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
