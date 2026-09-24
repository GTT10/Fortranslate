# Serial AMR LES Mach inlet integration — 2026-09-25

Implementation `c1b5e1f`, based on `363936f`.

## Change

The common/subcycle timestep selectors and common/per-level LES steppers accept
optional Mach and inlet-pressure arrays. They forward those controls and inlet
temperature to the shared viscous ghost sampler, including Dynamic's wider
filter halo. Recursive stepping forwards them through every interval. Ordinary
reservoir faces retain their packed-state requirement. Existing wall treatment,
coarse interpolation, integrated fluxes, reflux and transactional publication
remain unchanged. MPI LES and coupled application integration remain pending.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug. Fixed unit target and selected H2O2 spray-tree
runtime built successfully after the line-wrap correction below. Selected
execution was not run. `spray_les_amr_subcycle` passed: **1/1 tests**, no failures
or skips, 159.92 s (total 159.96 s).

All four models (Smagorinsky, WALE, Vreman and Dynamic) exercise new Mach 0.2,
700 K inlet cases with extrapolated or 120000 Pa pressure and moving isothermal
walls. The recursive timestep selector accepts both modes. Root stepping matches
the common-step oracle at existing tolerances. Nested ratios two and three retain
1/4/36 counts, composite conservation, synchronized parent cells, positive
temperatures, species closure and large-step rollback of state and diagnostics.
Existing periodic/wall/outlet cases remain included. No tolerance or timeout
changed. Provenance, Python syntax, project contract and whitespace passed.
This is focused numerical evidence, not turbulence or nozzle validation.

## Commands and artifacts

Logs: `build/mach-inlet-amr-les-checks/`. The initial `build-fixed.log` records
a line-length compile error; the continued argument line was wrapped before
`build-fixed-final.log`. Selected build output is in `build-selected.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_amr_subcycle --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_amr_subcycle$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Environment overrides are host-local workarounds, not application defaults.
