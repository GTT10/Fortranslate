# Two-level AMR default viscosity — 2026-09-25

Code `21888a0`, parent `352cf7b`. Serial two-level hydro resolves omitted
`difmag` to 0.1 and passes the resolved value to coarse/fine updates. MPI sparse
hierarchy, physical/periodic coarse stepping and private fine stepping likewise
resolve 0.1. Explicit zero remains available. The fine-step forwarding is
necessary because its PCM/PLM face-flux helper only applies an explicitly
supplied coefficient. Defaults are independent of solver/reconstruction.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Python 3.12.3, Ninja Debug; reused
`build/spray-amr-reactive-debug`. Final eight-rank CTest **1/1 passed**, 10.72 s
(10.78 s total), no skips. An initial run failed serial/MPI periodic parity:
the MPI fine helper had not received the resolved default. Forwarding the value
fixed the mismatch without changing tolerances. Both logs are retained.
Compiler and nonfatal display-authorization warnings remain.

PCM/PLM periodic fixtures compare omitted MPI against explicit serial 0.1;
physical interior/wall-touching fixtures reverse the assignment. State,
temperature and reflux comparisons match serial exactly; existing physical
boundary budget bound 2e-11 and failure checks pass. Separate reacting/full,
PPM, selected and full-matrix tests were not rerun. No external full-solver
parity or physical validation is claimed. Provenance, Python syntax, contract
and whitespace checks passed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_reactive_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^mpi_amr_default_hierarchy_3d_np8$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/two-level-default-difmag-checks/` contains `build.log`,
`rebuild.log`, `fix-build.log`, `LastTest-first.log` and `LastTest.log`.
