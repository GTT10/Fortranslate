# Pinned CTU physical boundaries — 2026-09-26

Implementation `bc7b5a29`, parent `3ed36cea`, local main, version 0.246.0.

The compiled reference-body oracle now includes unchanged pinned `PLM.H`
slopes/tracing and `pc_low_order_boundary`. It compares the complete zero-source
twelve-solve pipeline and final viscosity/species normalization on 4x4x4 grids.
Eight additional cases use independent mixture variation, dt=1e-5, flattening
enabled, and viscosity 0/0.1 for outflow, slip wall, no-slip wall and inflow.
Each case uses the same boundary type on all six faces; inflow uses one constant
reservoir state. Ghost construction is independently implemented in the driver.

The first sixteen cases passed. Adding inflow exposed a case-16 component-scaled
error of **0.003198352347347108**: Fortran normalized independently reconstructed
PLM species weights before the Riemann solve. Regular PPM inflow replacement now
retains these weights, matching the reference. Other PLM callers retain their
existing default normalization. No tolerance or golden data changed.

## Results

After the repair, **18 cases / 51,840 final-flux values passed**, maximum
component-scaled error **2.988913928686505e-15**, unchanged bound `3e-12`.
Related normal/halo units passed **2/2**, respectively 6.38 and 6.68 s,
13.11 s total; no skips. Provenance, Python syntax, project contract and
whitespace checks passed. GNU Fortran/C++ 13.3.0, Python 3.12.3, existing Ninja
Debug tests-enabled build. Existing compiler warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_normal_3d test_reactive_ppm_halo_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tools/check_pinned_ppm_ctu.py --reference /mnt/sn850/Projects/gtt2/external/PeleC --executable build/spray-amr-reactive-debug/tests/test_reactive_ppm_normal_3d --work-dir build/ppm-pinned-boundary-checks/repaired-oracle --shock-viscosity --physical-boundaries
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(unit_reactive_ppm_normal_3d|unit_reactive_ppm_halo_3d)$' --output-on-failure
python3 -m py_compile tools/check_pinned_ppm_ctu.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/ppm-pinned-boundary-checks/`, including `build-repair.log`,
`comparison-repaired.log`, `ctest.log`, `unit.log`, `provenance.log` and
`repaired-oracle/` sources, compiler log, input/output tables, hashes and report.
The earlier successful sixteen-case `oracle/` and failed `inflow-oracle/` outputs
are retained, with `comparison.log` and `comparison-inflow.log`.
Reference commit: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

This compares compiled reference bodies with a NASA7 SI adapter. It does not
establish full AMReX/PeleC executable parity. WENO, forcing, mixed boundary
types, pressure-outlet extensions, EB, selected mechanisms, MPI/application
execution, installed audit and physical validation remain outside this check.
