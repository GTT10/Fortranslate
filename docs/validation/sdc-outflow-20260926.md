# SDC outflow library route — 2026-09-26

Implementation `3d2b6cd9`, parent `90f56a2c`, local main, version 0.246.0.

`advance_reactive_sdc_3d` aliases the existing periodic API and accepts optional
boundary faces and a boundary-change inventory. Periodic and zero-gradient
outflow faces can be mixed. The outflow path fills the old state halos with
the existing boundary routine, extrapolates source rates at outflow faces,
wraps periodic axes and calls the existing halo CTU with outflow Riemann flags.
The final iteration's face fluxes determine the boundary inventory, separately
from volume-source changes. Existing source/reaction iteration is shared.

Molecular transport currently uses a periodic flux API, so a physical-boundary
call with built-in transport is rejected. Inflow, walls and prescribed outlet
pressure are not exposed by this SDC route. Application input remains periodic
until its boundary dispatch is connected.

## Completed checks

- GNU Fortran 13.3.0 Debug, fixed tests-enabled build
  `build/spray-amr-reactive-debug`, target `test_reactive_ppm_step_3d` passed.
- `unit_reactive_ppm_step_3d`: **1/1 passed**, 21.19 s test / 21.23 s total.
  Added nonuniform 4³ cases with x-outflow and periodic y/z. Zero-source SDC
  matches existing boundary CTU within 1e-13 scaled state error. Boundary
  inventories agree within 1e-10 absolute. With external force, all component
  inventories close after subtracting boundary and volume sources within
  1e-12 scaled error. Existing periodic reacting and source checks also pass.
- Provenance, Python syntax, project contract and whitespace passed. No failed
  or skipped tests; existing build warnings remain. The final post-build code
  edit only updated the routine's descriptive comment.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_step_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_step_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/sdc-outflow-checks/` contains build, provenance, CTest and
LastTest logs. This is library numerical evidence, not an external reacting
solver comparison or physical boundary-flow validation. Reacting outflow,
transport boundaries, application routing and AMR/MPI SDC remain open.
