# EB library reference viscosity default — 2026-09-25

Code `a33c631`, parent `499808a`. EB face fluxes, serial two-level/multipatch/tree
hydro, MPI common/recursive/full tree updates and MPI EB spray coupling resolve
omitted `difmag` to 0.1. Full/Strang and serial spray wrappers forward omission
into these paths. Shared private transport operators retain zero when invoked
for molecular transport alone. Explicit zero remains available.

An external EB halo now needs complete edge/corner state and matching volume
fractions when the coefficient is omitted, as already required for positive
explicit values. Face-only callers must supply explicit zero or upgrade their
halos. Integrated hierarchy callers request complete halos automatically.
Default changes do not resolve the PCM ancestor-geometry approximation.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Python 3.12.3, Ninja Debug; reused
`build/spray-amr-reactive-debug`. **3/3 CTests passed**: coupled two-level
17.65 s, coupled multipatch 12.32 s, three-rank MPI common tree 1.83 s
(31.85 s total). No failures or skips. Compiler warnings remain. Provenance,
Python syntax, contract and whitespace checks passed.

Coupled serial fixtures now leave their optional coefficient unallocated by
default and exercise omitted-argument common/recursive hydro, chemistry and
transport, conservation, operator-disabled parity and rollback. MPI common
hydro omits the coefficient while its serial reference explicitly supplies 0.1;
states/temperatures/boundary ledgers match exactly under existing layouts and
ownership modes. Existing explicit invalid/mismatched coefficient checks run.

MPI recursive/full and particle-coupled omitted paths compiled but were not
independently run here. Historical face-only callers, full regression matrix,
external geometry parity and physical validation were not audited. No tolerances
or golden hashes changed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_hydro_3d test_eb_two_level_full_3d test_eb_multipatch_full_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(unit_mpi_eb_tree_default_3d_3|unit_eb_two_level_full_3d|unit_eb_multipatch_full_3d)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/eb-library-default-checks/` contains `build.log`, `rebuild.log`
and `LastTest.log`.

The subsequent [face-only caller audit](eb-face-halo-default-20260925.md)
updates the two legacy halo fixtures and records three focused passing tests.
