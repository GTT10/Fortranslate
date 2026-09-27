# Embedded-wall heating application controls — 2026-09-27

Implementation: `96885be4`. The MPI EB tree app exposes
`embedded_wall_temperature` and forwards it through gas-only and particle common/
recursive SDC paths. Positive values require SDC, transport and conduction;
zero preserves the previous behavior. The serial app rejects enabling it.
The existing six-layer geometry cache also supplies the two-layer heat cache
after layout changes. Each interval is capped by the thermal frequency even
with automatic stepping disabled. Checkpoints bind the enabled wall temperature
through `PELEF_EB_TREE_EMBEDDED_HEAT_V1`; disabled checkpoints retain their old
header sequence. Selected-runtime targets now include all three new wall modules.

The checked-in [case](../../cases/spray_eb_wall_heat_3d/README.md) is a short
heated plane with heptane parcels. It does not claim full pinned Spray-EB parity,
embedded no-slip viscous traction or physical validation.

## Build and source identity

GNU Fortran 13.3.0 / OpenMPI 4.1.6, existing Debug configuration
`build/heptane-simple-mpi-debug`, selected pinned heptane Simple transport,
spray/MPI ON, tests/SUNDIALS OFF. The following targets built successfully:

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray pelef_mpi_reactive_eb_tree_3d_selected pelef_reactive_eb_tree_3d_selected --parallel 2 > build/eb-wall-app-build3.log 2>&1
python3 tools/check_mechanism_provenance.py --bundle build/heptane-3sp-simple-reference/heptane_3sp.json
python3 -m py_compile build/eb-wall-app/run_checks.py build/eb-wall-app/finish_checks.py build/eb-wall-app/check_defaults.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Spray executable SHA-256:
`5a8d065220fc80c038bdaa95d3d594ef0d74508bbbe8860a147bb8c7b2ca9b0a`.
The first build exposed omitted selected-runtime sources; the second exposed
two overlong wrapper declarations. Both were corrected; failure logs
`build/eb-wall-app-build.log` and `build/eb-wall-app-build2.log` are retained.

## Bounded application checks

Scripts and all inputs/logs are under `build/eb-wall-app/`:
`run_checks.py`, `finish_checks.py`, `check_defaults.py`; reports `result.json`
and `defaults-result.json`. Runs use `/usr/bin/mpiexec -n N` with DISPLAY,
LD_LIBRARY_PATH and LIBRARY_PATH removed, HWLOC_COMPONENTS=-x86 and
OMPI_MCA_hwloc_base_binding_policy=none. Each case directory contains `tree.nml`,
`spray.nml`, `level_set.dat` and `fortran.log`.

Completed checks:

- `continuous`: one rank, 4³ cells, two 100 ns SDC steps, 450 K embedded wall,
  initially 400 K gas, four parcels. Wall energy `9.561100114980793e-9 J`;
  total-energy budget residual `-1.897353801849633e-19 J`.
- `checkpoint` / `restart`: stop after step one, then resume with two ranks.
  Final gas CSV bytes and sorted per-rank particle rows match the continuous run.
- `changed`: changing wall temperature to 451 K is rejected on two-rank restart
  before final CSV publication.
- `cap2`: with automatic stepping disabled and thermal CFL `1e-8`, proposed
  `1e-7 s` is capped to `1.300563389122448e-11 s`. Checkpoint stop after one step
  exits successfully. The earlier `cap` attempt intentionally used only one
  allowed step but lacked checkpoint-stop; it hit the normal maximum-step error.
  It is retained and not counted as a completed run.
- `default`: disabled new setting reproduces `build/spray-face-transport`
  gas/particle CSV bytes exactly.
- `recursive`: the heated single-level particle application completes with
  recursive stepping enabled. Actual multilevel recursive gas conservation is
  covered by the preceding [integrator test](eb-wall-sdc-integrator-20260927.md),
  not by this single-level application case.

No tolerance changed. Expected compiler/linker/runtime warnings remain in logs;
the corrected positive runs exit zero. No full CTest matrix, installed-binary
audit or full reference-solver/physical comparison was performed. The original
full positive-rate decane jet executable was untouched and remained running;
its live state reached 0.3458 ms / 1,922 steps during this work.
