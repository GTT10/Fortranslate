# Embedded-wall heat source in SDC transport — 2026-09-27

Implementation: `9dfb15c9`. This connects the
[quadratic stencil](eb-boundary-gradient-20260927.md) to conserved-state
temperature recovery and the existing molecular transport database, including
its Simple-fit dispatch. Conductivity is evaluated at the cut cell, as in the
pinned scalar boundary-flux application. The geometry's integrated normal gives
the reference aperture-vector area; physical boundary centroids are converted
to cell-relative coordinates.

`compute_eb_wall_heat_source_3d` accepts a prepared geometry/state containing
donor cells and an interior index box. It returns energy per fluid volume per
second and the signed total heat added to gas. Covered or missing active donors
reject the complete operation; outputs remain unchanged. The returned frequency
is an absolute row-sum bound for frozen conductivity divided by fluid-volume
heat capacity (`rho*cv`). It is not a nonlinear positivity guarantee.

`compute_sdc_eb_transport_3d` now optionally adds this source to its energy rate
and returns time-integrated wall energy plus the frequency. Wall temperature,
wall energy and frequency arguments must be supplied together; conduction must
be enabled. Rate, face fluxes, minimum theta and wall outputs publish together
after all checks. The default call retains its previous behavior.

The current SDC entry requires active stencil donors inside its supplied patch.
It fails if donors are missing; it does not fabricate or silently omit them.
The prepared source API can consume a wider same-resolution geometry, but AMR
stage halo exchange, application controls, global wall-ledger accumulation,
timestep selection and restart identity are not yet wired to this new option.
This is library integration evidence, not a completed Spray-EB case.

## Completed checks

GNU Fortran 13.3.0; existing Debug core configuration
`build/heptane-simple-mpi-debug` (OpenMPI 4.1.6, spray ON, tests/SUNDIALS OFF).
No compiler/backend setting changed and the running decane executable was not
rebuilt. Core build log: `build/eb-wall-heat-core-build.log`.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_core --parallel 2 > build/eb-wall-heat-core-build.log 2>&1
/usr/bin/gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -I build/heptane-simple-mpi-debug/modules tests/unit/test_eb_wall_heat_source_3d.F90 build/heptane-simple-mpi-debug/libpelef_core.a -o build/eb-wall-heat-unit/core-check
build/eb-wall-heat-unit/core-check
/usr/bin/gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -I build/heptane-simple-mpi-debug/modules tests/unit/test_reactive_eb_transport_3d.F90 build/heptane-simple-mpi-debug/libpelef_core.a -o build/eb-wall-heat-unit/existing-transport-check
build/eb-wall-heat-unit/existing-transport-check
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

All commands passed. Two focused unit executables passed:

- New source test: 6³ plane geometry, quadratic nitrogen temperature field.
  Analytic wall heat `-0.11456198911385618 W`; integrated source
  `-0.11456198911385673 W`. The test checks source support, positive frequency,
  exact SDC wall-energy routing, missing-ledger rejection and a late invalid
  donor with full output rollback.
- Existing 3D EB transport test: axis cases, full reactive driver and curved
  geometry cases passed with their existing tolerances. Reported axis total
  conservation errors range from `7.1e-16` to `1.4e-15`.

An earlier direct compile of the changed sources also passed the new test;
it is not counted as another distinct test. CTest registration was added, but
no full CTest matrix or installed-binary/reference-solver audit was run.
No Python source changed in this increment.
