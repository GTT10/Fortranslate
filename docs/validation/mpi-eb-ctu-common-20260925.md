# MPI EB tree CTU common intervals — 2026-09-25

Implementation `2cb5317`, parent `0f7df17`; pinned PeleC reference
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8` is unchanged.

The common-interval MPI hydro entry now accepts the serial six-layer extended
geometry cache. Ranks collectively agree on its presence and all geometry
metrics, then verify physical centers. Only owners construct wide metric fields
and advance CTU patches. Existing distributed peer reconciliation, reflux,
restriction and boundary/wall accounting publish the accepted result together.
The geometry validator accepts collectively agreed padding; parent/child physical
consistency remains checked by the ordinary physical-geometry validation.

## Focused checks

GNU Fortran 13.3.0, OpenMPI through `/usr/bin/mpiexec` and `/usr/bin/mpifort`,
Python 3.12.3, Ninja Debug; existing fixed-mechanism build
`build/spray-amr-reactive-debug`. No compiler or MPI configuration changed.

- `unit_mpi_eb_tree_hydro_3d_3`: PASS, 1.75 seconds.
- `unit_mpi_eb_tree_ctu_hydro_3d_3`: PASS, 6.49 seconds.
- Total: 2/2 passed in 8.29 seconds, no skips or interruptions.

CTU checks one five-node, three-level z-plane tree at 1000 K with discontinuous
densities and `difmag=.037`, under distributed ownership and all nodes owned by
the last rank (two idle ranks). State, temperature and boundary ledger match
serial exactly. Wall and composite conservation retain the `3e-14` threshold;
maximum MPI CTU composite residual is `4.634e-15`. Initial velocities are zero;
this does not isolate a nonzero viscosity correction.

Rejection tests cover solver, timestep, bad density, viscosity mismatch,
physical halo context, cache-presence mismatch and a rank-local extended face
centroid change. State/temperature and output ledgers remain unchanged.

The initial build failed on a new overlong test line, corrected before testing.
Initial CTU testing passed, but legacy zero-viscosity hydro failed: its wrapper
omitted `difmag` in the zero branch and inherited the changed nonzero default.
The wrapper now forwards the explicit coefficient in both branches. Both tests
passed after this correction. Compiler warnings and X authorization messages
from the legacy MPI test environment remain; no warning-clean claim is made.

Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_hydro_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_(ctu_hydro_3d_3|hydro_3d_3)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, contract and whitespace checks passed. Mechanism SHA remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Artifacts in `build/mpi-eb-ctu-checks/`: `build.log`, `build-final.log`,
`build-zero-fix.log`, `LastTest-initial.log`, `LastTest-final.log`.

MPI recursive CTU, coupled MPI dispatch, tree application selection/restart,
other physical boundaries, external solver parity and physical validation remain
separate work. This is owner-local numerical integration evidence, not whole
PeleC qualification.
