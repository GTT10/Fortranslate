# MPI spray library default viscosity — 2026-09-25

Code `5cc224a`, parent `9a0d583`. MPI regular spray hydro and mixture coupling
resolve omitted `difmag` to 0.1, independent of solver/reconstruction. The
mixture wrapper forwards that value to gas hydro. Explicit zero and existing
collective validation remain. This aligns library calls with application inputs.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Python 3.12.3, Ninja Debug; reused
`build/spray-amr-reactive-debug`. Final **2/2 four-rank CTests passed**: PPM
9.84 s, PLM 15.17 s (25.07 s elapsed), no failures or skips. Compiler warnings
remain. Provenance, Python syntax, contract and whitespace checks passed.

Each test runs four policies and eight boundary configurations on two-cell
slabs. The new HLLC policy omits the MPI coefficient using an unallocated
optional scalar and compares against explicit serial 0.1. Conserved state and
temperature agree exactly; boundary budgets pass the existing 2e-13 scaled
bound. Explicit zero/0.1/0.037 policies and collective CFL, selector and
coefficient failure/rollback checks also pass.

An earlier passing run still supplied an explicit coefficient in the added
policy. Inspection caught this test omission; calls were corrected and both
tests rerun. That earlier log is not counted as omitted-default evidence.
The particle-coupled omitted-default path compiled but was not independently
run; these are hydro library checks, not physical validation or a full matrix.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_ppm --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_(ppm|plm_viscosity)_4$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/mpi-spray-library-default-checks/` contains `build.log`,
`rebuild.log`, `LastTest-explicit.log` and final `LastTest.log`.
