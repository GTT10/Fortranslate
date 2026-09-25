# Recursive and coupled MPI EB CTU — 2026-09-25

Implementation `8bc58ac`, parent `702a0da`; PeleC pin remains
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

Recursive MPI hydro now consumes the optional extended geometry cache, with
collective presence/metric validation and physical-center checks. Wide halos
use ancestor time interpolation at each substep. Owners advance CTU patches;
existing distributed synchronization accumulates interval fluxes, refluxes and
restricts. The coupled MPI entry forwards the cache to both common and recursive
hydro inside transactional R/2 T/2 H T/2 R/2 composition.
The recursive zero-viscosity PCM branch now forwards its explicit coefficient
instead of inheriting the lower-level nonzero default.

## Completed checks

GNU Fortran 13.3.0, OpenMPI `/usr/bin/mpiexec` and `/usr/bin/mpifort`, Python
3.12.3, Ninja Debug, existing `build/spray-amr-reactive-debug` configuration.

| CTest | Result | Seconds |
| --- | --- | --- |
| `unit_mpi_eb_tree_subcycle_hydro_3d_3` | PASS | 4.20 |
| `unit_mpi_eb_tree_ctu_subcycle_hydro_3d_3` | PASS | 17.58 |
| `unit_mpi_eb_tree_ctu_full_3d_3` | PASS | 106.76 |
| `unit_mpi_eb_tree_ctu-common_full_3d_3` | PASS | 84.15 |

4/4 passed, total 212.72 seconds; no failures, skips or interrupted runs.
CTU uses the five-node, three-level z-plane density-discontinuity fixture with
1000 K fluid temperatures and `difmag=.037`. Recursive hydro checks distributed
and all-on-last-rank ownership, including two idle ranks. Coupled cases check
explicit chemistry with distributed ownership and implicit chemistry with all
nodes on the last rank. Both transport and chemistry are enabled; existing
non-CTU switch coverage was not rerun as a full matrix.

State, temperature, boundary ledger and limiter match serial exactly. Chemical
and wall ledgers and composite inventory retain the `3e-14` threshold. Maximum
MPI CTU residuals: hydro `4.468e-15`, coupled recursive `4.303e-15`, coupled common
`4.303e-15`. Initial velocities are zero; these are not isolated nonzero-viscosity
oracle checks. Recursive hydro rejects inconsistent dt, coefficient, geometry
context/cache and bad density collectively, preserving inputs/outputs. Coupled
late hydro rejection preserves state and diagnostics after private earlier stages.

Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_subcycle_hydro_3d test_mpi_eb_tree_subcycle_full_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_(subcycle_hydro_3d_3|ctu_subcycle_hydro_3d_3|ctu(-common)?_full_3d_3)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, Python syntax, project contract and whitespace checks passed.
Mechanism SHA unchanged:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Artifacts: `build/mpi-eb-ctu-recursive-checks/build.log` and `LastTest.log`.
Compiler warnings remain. Tree application selection/restart, broader physical
boundaries, whole-solver reference parity and physical validation remain separate.
