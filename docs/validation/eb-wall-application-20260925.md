# Reflected EB wall geometry and thermal application boundaries — 2026-09-25

Implementation `c854c30`, based on `7a951c9` on `main`.
PeleC pin remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

Configured plane, extruded and nodal geometry now reflects outside fixed domain
walls and extrapolates at other faces. Nodal donor mapping supports repeated
reflection across a domain narrower than the halo and mixed wall/outflow faces.
The physical geometry remains unchanged. A reflected plane is reconstructed
from nodal values when its normal direction has a reflecting wall.

Fixed and selected single-level gas applications accept `slip_wall`,
`no_slip_wall` and `isothermal_wall`, plus six wall temperatures. Isothermal input
requires transport and conduction. The application invokes the aperture/volume
weighted wall thermal timestep bound and forwards temperatures into the full gas
step. Checkpoint identity adds `wall-v1` and temperatures for isothermal cases;
legacy and inflow-only identity text remains unchanged. The shared reader opt-in
still prevents unsupported spray/tree front ends from silently accepting walls.
An example is `cases/reactive_eb_3d/plm_ctu_wall.nml`.

These controls concern fixed Cartesian domain walls, not thermal conditions on
the embedded cut surface. Moving walls and periodic transport remain pending.

## Focused evidence

GNU Fortran 13.3.0, Python 3.12.3, Ninja Debug, existing OpenMPI configurations.
Fixed tests-on: `build/spray-amr-reactive-debug`. Selected tests-off:
`build/spray-amr-reactive-selected-debug`, `mechanisms/h2o2_full.json`.

- Final unit run **3/3 passed**, 15.32 s total: config 0.40 s, checkpoint 0.18 s,
  transport halo 14.69 s. No failures or skips.
- Geometry checks cover legacy extrapolation, six reflecting faces and mixed
  lower-wall/upper-outflow faces across plane/extruded axes and a 3D nodal field.
  They check the physical crop, repeated state/metric donor correspondence and
  edge/corner nodal values. The coupled wall suite now uses the configured
  geometry builder rather than constructing its own reflected fixture.
- Fixed and selected application checkers each completed **5 successful runs and
  6 expected rejections**. Inert/reacting adiabatic cases and coupled isothermal
  cases pass inventory checks. Continuous and checkpoint-resumed final CSVs
  match byte-for-byte. Wall-temperature changes reject at restart without
  creating the requested CSV. Invalid temperature, disabled conduction and
  incompatible methods reject during input validation.

The wall-rate formula was isolated in the preceding timestep unit check; these
short application runs do not separately demonstrate a thermal-dominated step.
The geometry fixtures are planar fields through the three input representations,
not arbitrary curved-wall physical validation. No tolerances were relaxed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_eb_3d test_reactive_eb_3d_config test_eb_transport_halo_3d test_reactive_eb_3d_checkpoint --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_(reactive_eb_3d_(config|checkpoint)|eb_transport_halo_3d)$' --output-on-failure
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 python3 tools/check_eb_ctu_application.py --executable build/spray-amr-reactive-debug/pelef_reactive_eb_3d --work-dir build/eb-wall-application-checks/fixed --wall
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 python3 tools/check_eb_ctu_application.py --executable build/spray-amr-reactive-selected-debug/pelef_reactive_eb_3d_selected --work-dir build/eb-wall-application-checks/selected --wall --selected
python3 -m py_compile tools/check_eb_ctu_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, contract and whitespace checks passed. Compiler warnings
remain, including the existing internal test-provider executable-stack warning.
No build or runtime failures occurred. Mechanism SHA remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Artifacts: `build/eb-wall-application-checks/{build.log,build-selected.log,LastTest.log}`,
fixed `fixed/run-z2irk_pg`, selected `selected/run-g8qfumix`.
Earlier geometry-only evidence is in `build/eb-reflected-geometry-checks/`.
No full matrix, installed-binary audit, external-PeleC comparison or physical
wall-flow qualification was performed here.
