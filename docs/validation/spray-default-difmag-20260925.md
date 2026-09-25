# Regular spray reference viscosity default — 2026-09-25

Code `94dff7a`, parent `ebbab3b`. Serial and MPI regular spray applications
resolve omitted/-1 `difmag` to 0.1 for every supported reconstruction and solver,
matching the pinned PeleC parameter default. Explicit nonnegative overrides
remain available. Existing checkpoint contexts bind resolved values; old
zero-policy continuations must specify zero, while automatic and explicit 0.1
are equivalent. No golden hashes or tolerances changed.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Python 3.12.3, Ninja Debug; reused
`build/spray-amr-reactive-debug`. Fixed serial/MPI executables built.
**2/2 CTests passed**: serial PLM input/restart 27.62 s and MPI PLM input/restart
9.29 s (36.95 s total), no failures or skips. Compiler warnings remain.
Provenance, Python syntax, project contract and whitespace checks passed.

The serial Mach-inlet coupled fixture and MPI periodic-mixture fixture verify
automatic versus explicit 0.1 restart equivalence. MPI continuation changes
2 to 4 ranks. Existing positive coefficient, zero, changed-coefficient,
invalid-input and particle/gas continuation checks pass. This increment does
not independently rerun selected variants, serial PCM, PPM, other solvers or
the full suite. This is bounded restart/application evidence, not physical
validation. Direct-library and EB defaults still need unification.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_spray_3d pelef_mpi_spray_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_(mpi_)?plm_difmag_restart$' --output-on-failure
python3 -m py_compile tools/check_spray_boundaries.py tests/spray/test_mpi_spray_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/spray-default-difmag-checks/build.log`, `LastTest.log`.
