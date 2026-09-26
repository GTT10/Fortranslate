# Transactional drift after particle halo kicks — 2026-09-27

Implementation `c38f9d92`, based on `ce7c8c5b`, local `main`, version 0.246.0.
Pinned PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Implementation

`kick_mpi_amr_mixture_halo` accepts optional `drift_interval`, `boundary_faces`,
`wall_change` and `boundary_change` as one group. It first computes the existing
private kick, then drifts with the resulting velocity. Durations are independent.
The existing level-set drift handles the supplied EB geometry or an all-fluid
regular domain, with a limit of 32 contacts/crossings per particle. Physical
boundary labels must agree with sampling periodicity and across ranks.

EB impulse/work and physical boundary inventories are separately reduced over
MPI ranks. Outflow also exports each component's vapor-species inventory and
zeros component masses. Dead entries remain in the population. Source deposition
uses the pre-drift sampling positions. Source, particles, components and both
ledgers publish together only after all operations succeed. IDs and assigned
nodes are retained; ownership migration, ghost regeneration, breakup, population
accounting and recursive application scheduling remain caller work.

## Completed verification

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6. Final `spray_mpi_amr_step_2` and `_4` passed
(2/2, 0.81 s and 0.78 s; total 1.64 s).

The extended fine-halo fixture checks a zero-kick periodic drift of 2e-4 s
against wrapped ballistic coordinates within 1e-15, with unchanged velocity
and components and zero external exchange. A one-second drift exceeds the
crossing budget after a private nonzero kick and rejects; checked particle
positions/velocities, components, source arrays and sentinel ledgers remain
unchanged. A regular outflow case zeros parcel/component mass and exports
total mass and each vapor-species inventory exactly. Existing halo-kick,
composite coupling and capped-view assertions run in the same executables.

Provenance preflight, build, Python syntax, project contract and whitespace
checks passed. No build/test failures or skips occurred. An earlier two-test
run passed before adding the outflow assertions; the timings above describe
the final run only. EB contact and reflecting-boundary paths were compiled
through the existing drift routine but were not independently exercised here.
No full matrix, installed audit, external solver parity, recursive application
or physical validation was performed.

## Reproduction

Artifacts: `build/spray-halo-drift-checks/`, including final CTest details.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_step --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_step_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
