# Periodic EB geometry and coupled application input — 2026-09-25

Implementation `e9782c8`, based on `2b9adde` on `main`.
PeleC pin remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

Configured geometry now wraps outside paired periodic faces for plane, extruded
and nodal inputs. A shared seam check compares opposing apertures and tangential
face centroids. With an extended geometry it additionally checks wrapped halo
cell volume fractions and centroid offsets against physical donors. Configuration
construction publishes only compatible geometry; full/Strang CTU drivers also
check the extended metric seams before stepping. Transport uses the shared seam
check instead of its earlier aperture-only check.

The full gas driver now accepts periodic built-in transport donors. Fixed and
selected single-level applications accept paired periodic labels and use the
existing boundary-bound checkpoint identity. Unpaired input, incompatible periodic
geometry and changed periodic restart settings reject. Custom transport providers
remain limited to inflow/outflow. An example is
`cases/reactive_eb_3d/plm_ctu_periodic.nml`.

## Focused evidence

GNU Fortran 13.3.0, Python 3.12.3, Ninja Debug, existing OpenMPI configurations.
Fixed tests-on: `build/spray-amr-reactive-debug`. Selected tests-off:
`build/spray-amr-reactive-selected-debug`, `mechanisms/h2o2_full.json`.

- **3/3 unit tests passed**, 16.94 s total: configuration 0.53 s, checkpoint
  0.14 s, transport halo 16.22 s. No failures or skips.
- Geometry fixtures include repeating extruded and 3D nodal cut fields, with
  six halo cells around a four-cell physical width. Repeated wrapping at edges
  and corners preserves donor volume fractions. A nonperiodic plane geometry
  rejects without replacing the prior output.
- Full-fluid and tangentially periodic cut-cell library cases now advance
  reaction, transport and CTU with inventory closure at the existing `2e-12`
  scaled criterion. Existing periodic species-limiter and regular-flux checks
  also pass.
- Fixed and selected application checkers each completed **5 successful runs and
  6 expected rejections**. Inert, reacting, coupled, checkpoint-stop and resume
  runs pass inventory checks. Continuous/resumed final CSVs match byte-for-byte.
  Changed periodic pairs reject at restart without creating the requested CSV;
  unpaired labels reject in configuration validation; a plane cut normal to the
  requested periodic axis rejects during geometry construction.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_eb_3d test_reactive_eb_3d_config test_eb_transport_halo_3d test_reactive_eb_3d_checkpoint --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_(reactive_eb_3d_(config|checkpoint)|eb_transport_halo_3d)$' --output-on-failure
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 python3 tools/check_eb_ctu_application.py --executable build/spray-amr-reactive-debug/pelef_reactive_eb_3d --work-dir build/eb-periodic-application-checks/fixed --periodic
env -u DISPLAY -u LIBRARY_PATH -u LD_LIBRARY_PATH HWLOC_COMPONENTS=-x86 python3 tools/check_eb_ctu_application.py --executable build/spray-amr-reactive-selected-debug/pelef_reactive_eb_3d_selected --work-dir build/eb-periodic-application-checks/selected --periodic --selected
python3 -m py_compile tools/check_eb_ctu_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, contract and whitespace checks passed. No build/runtime
failures occurred. Compiler warnings remain, including the internal test-provider
executable-stack warning. Mechanism SHA remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Artifacts under `build/eb-periodic-application-checks/`: `build*.log`,
`LastTest.log`, fixed `fixed/run-0r8n3vu2`, selected `selected/run-byz016l0`.
No full matrix, installed-binary audit, external-PeleC comparison or physical
validation was performed. Independent inlet mixtures and spray/tree/MPI boundary
dispatch remain pending; this is not complete PeleC parity.
