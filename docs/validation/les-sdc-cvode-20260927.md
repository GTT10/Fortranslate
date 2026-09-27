# Regular LES SDC forced-reactor dispatch, 2026-09-27

## Source and scope

Implementation: `c11e27b83aac746f53609b8ea882caf3d62f83d6`, based on
`45e126ffa1ecdaf0c492a9cb0e7080d5ccbcb169` on local `main`.
`advance_les_sdc_3d` now forwards an optional `forced_reactor_integrator` to
the regular gas SDC driver. Omitting it retains the native reactor. The
existing candidate-state, reaction-history, source-estimate, inventory and
face-flux transaction surrounds callback execution.

This is library integration. Regular gas application LES input routing remains
open. AMR and MPI particle applications use their separately connected paths.
No external PeleC parity, physical validation, selected-mechanism matrix,
installation audit or changed-rank run is claimed here.

## Completed checks

Reused `build/amr-cvode-app-debug`: GNU Fortran 13.3.0, Debug, MPI/spray/SUNDIALS
enabled, Open MPI 4.1.6, SUNDIALS 7.2.0 static Fortran interfaces. Compiler and
MPI paths are `/usr/bin/gfortran` and `/usr/bin/mpifort`; SUNDIALS package path
is `build/sundials-7.2-fortran/install/lib/cmake/sundials`. No compiler or
backend configuration changed. These tests execute the fixed-mechanism serial
library. Existing real-comparison and executable-stack warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/amr-cvode-app-debug \
  --target test_les_sdc_cvode_3d test_les_sdc_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 \
  OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/amr-cvode-app-debug \
  -R '^les_sdc(_cvode)?_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, build, Python syntax, project contract and whitespace checks passed.
The pinned YAML SHA-256 is
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
CTest: **2/2 passed, 0 failed, 0 skipped**, total 13.48 s:

- `les_sdc_3d`: 10.98 s; existing native/inert four-closure, boundary,
  flux reconstruction, supplied-halo and rollback checks.
- `les_sdc_cvode_3d`: 2.45 s; periodic 4³ cells, nonuniform velocity/pressure,
  H2/O2/N2 at 1100 K, molecular transport, forcing and two SDC iterations of
  a 1e-10 s interval. Smagorinsky, WALE, Vreman and Dynamic each invoke real
  CVODE 128 times, produce positive diffusivity and nonzero reaction history,
  and agree with the early-time native reference to a maximum normalized
  state difference of `2.0261570199409107e-15` (gate `2e-7`). Maximum normalized
  Euler inventory residual is `1.0821484339580831e-15` (gate `3e-12`).
  Deliberate callback failure on call 65 corrupts private arguments and confirms
  exact preservation of gas, temperature, reaction history, external source
  estimates, boundary/source inventories, diffusivity and all three face-flux
  arrays for each closure.

Artifacts: `build/les-sdc-cvode-checks/build.log`, `ctest.log`, `LastTest.log`.
No tolerances or golden data changed.
