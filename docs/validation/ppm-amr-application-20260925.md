# Coupled AMR gas application PPM — 2026-09-25

Source: base `1398d48`, implementation `144bcd7`. GNU Fortran 13.3.0,
Open MPI 4.1.6 (`/usr/bin/mpiexec`), Ninja Debug, Python 3.12.3. Existing fixed
MPI tests-enabled `build/spray-amr-reactive-debug` and selected H2/O2
tests-disabled `build/spray-amr-reactive-selected-debug` were reused separately.
Artifacts and exact inputs/checker: `build/ppm-amr-application-checks/`.

The fixed/selected serial/MPI patch-tree gas front ends now enable PPM in
input. The coupled split uses limited coarse donors for molecular transport
with PPM, matching the hydro choice. Mach inlet timestep evaluation uses
four-layer limited donors. Existing checkpoint context serializes the entire
flow configuration except continuation limits/output name, so reconstruction
selection is already bound to restart without a format change.

## Completed checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target pelef_amr_tree_reactive_3d pelef_mpi_amr_tree_reactive_3d --parallel 2`
  and `cmake --build build/spray-amr-reactive-selected-debug --target pelef_amr_tree_reactive_3d_selected pelef_mpi_amr_tree_reactive_3d_selected --parallel 2`: passed.
- `python3 build/ppm-amr-application-checks/check.py`: six completed output
  runs and one expected context rejection. Fixture derives from
  `cases/amr_tree_reactive_3d/hotspot.nml` and `branched.dat`, with full H2/O2
  or selected H2/O2, PPM and recursive subcycling. Chemistry and molecular
  transport are enabled with outflow and moving no-slip walls.
- Fixed serial versus two ranks, and selected serial versus two ranks:
  composite CSV and budget CSV are byte-identical at 1e-9 s (one root step).
- Fixed two-rank checkpoint continuation to 2e-9 s (two cumulative root steps)
  agrees byte-for-byte in composite/budget CSVs between serial and three ranks.
  This compares two continuations from the same checkpoint, not a separately
  rerun uninterrupted trajectory. Changing PPM to PLM rejects before overwriting
  a sentinel output file.
- Initial checker expected a final `Steps:` report after
  `stop_after_checkpoint=.true.`; the application correctly stopped before
  final output, causing a checker assertion. The test now writes the checkpoint
  without early stop. Initial logs are retained; the already completed fixed
  serial reference was reused. No solver/tolerance change addressed this.
- `python3 -m py_compile build/ppm-amr-application-checks/check.py tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, `git diff --check`: passed.

Compiler warnings remain. No full matrix, selected restart, automatic/scheduled
regridding, common-step application run, Mach inlet application run, convergence
or external/physical validation was rerun here. Library common/recursive hydro
and Mach-boundary checks are separately dated evidence. Exact reference
primitive-pressure parity remains open; this does not enable PPM in two-level
AMR or spray front ends.
