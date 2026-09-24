# EB Mach inlet timestep bounds — 2026-09-25

Implementation `f2bb4f9`, base `b55cee3`.

Serial and MPI EB tree timestep APIs now accept the inlet controls used by
their stepping APIs. When controls are supplied, current snapshot halos supply
additional signal-speed and molecular-diffusivity maxima. Only open face donors
are scanned; unused interior/edge/corner and covered-face halo padding is ignored.
The shared halo-bound helper evaluates transport coefficients through the
existing single-cell regular transport bound, then the existing EB metric/
volume bound applies these maxima to every advanced cell, including hidden
parents. Existing common and recursive root-interval scaling remains.

MPI constructs owner-local halos with collective inlet-control validation,
reduces donor maxima and chooses the global minimum interval. No gathered
global field is introduced. Invalid or rank-inconsistent input leaves the
caller's timestep and input fields unchanged. Without inlet arguments the
previous estimator path remains.

## Completed focused verification

GNU Fortran 13.3.0, Open MPI 4.1.6 with `/usr/bin/mpiexec`, existing Ninja Debug
fixed elementary H2/O2 build. Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_transport_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_transport_3d_2$' --output-on-failure --output-log build/eb-mach-timestep-checks/ctest.log
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 passed**, test **8.07 s**, total **8.11 s**. Existing x/y/z plane-cut,
four/five-patch transport conservation and serial/MPI checks remain. New
five-patch z-cut checks prescribe Mach 10, 2000 K and 101325 Pa at all faces.
For common/recursive stepping and hydro-only, mixed and diffusion-dominated
CFL settings, the inlet strictly tightens the bound and serial/MPI results are
exactly equal. Both round-robin ownership and all patches on the last rank are
tested. A one-rank 1 Pa pressure discrepancy rejects without changing the
supplied timestep. Existing field-preservation and inconsistent CFL/mode tests
remain. No tolerances changed.

Build, provenance, Python syntax, project contract and whitespace passed.
No failures or skips; compiler warnings remain. Artifacts:
`build/eb-mach-timestep-checks/{build.log,ctest.log,LastTest.log}`.
Application input/checkpoint wiring is still pending. No selected-mechanism
execution, full matrix, production scaling, external solver run or physical
nozzle validation was performed.
