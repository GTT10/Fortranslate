# EB isothermal wall timestep bound — 2026-09-25

Implementation `4a1e3e1`, based on `7ff0f17` on `main`.
PeleC pin remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

`reactive_eb_isothermal_timestep_3d` uses the same half-cell thermal law as the
wall flux. Each adjacent isothermal domain face contributes its regular-cell
thermal rate multiplied by aperture / fluid-volume fraction. Rates add at cell
edges and corners; the maximum cell rate sets `dt = cfl / rate`. Covered cells
and closed faces contribute nothing. No active thermal face returns an unlimited
bound. Invalid active wall data or arithmetic overflow rejects with `dt=0`.

This is a library timestep bound. Applications still need wall-compatible
extended geometry, wall input/restart binding and invocation of this bound.
It is not a change to the embedded cut-surface thermal model.

## Focused check

GNU Fortran 13.3.0, Python 3.12.3, existing Ninja Debug/OpenMPI configuration in
`build/spray-amr-reactive-debug`. `unit_eb_transport_halo_3d` passed **1/1** in
15.58 s (15.62 s CTest total), no failures or skips. It retains the coupled wall
flux/inventory checks and adds:

- Full-fluid and tangential-cut bounds matching the regular Fortran wall-corner
  rate at `3e-12` relative tolerance.
- A planar fluid layer occupying 20% of the last x cell, with a fully open heated
  upper x face: the allowed interval is 20% of the regular single-wall interval.
- Invalid temperature on a fully covered lower face is ignored; invalid active
  upper-face temperature rejects and resets the timestep to zero.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_transport_halo_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_transport_halo_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, Python syntax, project contract and whitespace checks passed.
Mechanism SHA remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Artifacts: `build/eb-wall-timestep-checks/{build.log,LastTest.log}`.
Compiler warnings remain, including the test's existing executable-stack warning.
No application, full matrix, installed-binary audit, external-PeleC comparison or
physical wall-flow validation was performed for this increment.
