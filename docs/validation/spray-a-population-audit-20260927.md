# Spray-A growing-population audit — 2026-09-27

The live-run status below is historical. The 12-event run was interrupted
after an input-normalization defect was found. Corrected completed runs and
child-release evidence supersede it in the
[normalization record](spray-a-khrt-threshold-20260927.md).

Implementation `14f7831e90b6f6cae0242d9184c7a7ef8fe4a53c`, based on
`6f56d1c3`. The MPI EB spray application now reports the accepted particle-ID
watermark and cumulative computational parcel splits. Both values already
persist in checkpoints; this change only exposes them in final output.
`check_jet_inventory.py --allow-khrt-children` requires KHRT preparation
metadata and those runtime records, verifies unique positive IDs within the
reported bound, and retains the existing fuel-balance gate (`1e-8`). It reports
IDs beyond injection and computational splitting separately; this is not a
count of surviving descendants or proof of breakup dynamics.

Spray-A preparation now computes maximum steps from the prepared duration,
fixed timestep and number of event boundaries. The rate table has nonuniform
interval widths, so an event-count-only bound is insufficient in general.
No physical coefficients, tolerances or numerical update formulas changed.

## Completed checks

GNU Fortran 13.3.0/Open MPI 4.1.6, existing clean Debug selected dodecane build:

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/dodecane-simple-mpi-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
python3 -m py_compile tools/check_jet_inventory.py tools/prepare_spray_a_reference.py build/run_spray_a_children.py
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Build, MPI ABI audit and listed static checks pass. Build log:
`build/spray-a-children-build.log`. The previous decane output passes the
unchanged default population bound and fuel gate
(`build/spray-a-children-decane-regression.json`). The new mode rejects old
Spray-A output lacking the required watermark, producing no result file
(`build/spray-a-population-rejection.json`). These checks do not establish
successful runtime child release.

## Extended run in progress

```sh
python3 tools/prepare_spray_a_reference.py --bundle build/dodecane-simple-reference/dodecane_lu.json --reference-source build/pelec-pinned-hit-source --output build/spray-a-children --events 12
python3 build/run_spray_a_children.py
```

This retains the previous physical setup, extends to twelve source rate
intervals (4.4 microseconds), and injects 48 parcels before any breakup-created
population. The one-rank driver and binary identity are recorded in
`build/spray-a-children/execution.json`. At this record, solver PID 676113 was
confirmed live (`Rl`, elapsed 5 min 20 s); the run is not counted as passed.
The execution handle is session 92434. Do not restart merely because a poll
times out. Upon successful completion, inspect the output population and run:

```sh
python3 tools/check_jet_inventory.py --run-dir build/spray-a-children --fuel-species NC12H26 --allow-khrt-children --output build/spray-a-children/inventory.json
```

Child release, the new runtime bound on a completed growing population,
full-duration reference parity and physical validation remain unproven here.
The separate frozen Jet continuation was also confirmed live as PID 648941;
its numerical binary was not changed by this work.

## Live progress observation

Subsequent read-only GDB attachment, followed by explicit detach, observed
accepted time 4.0 microseconds / step 22 / particle-ID watermark 1672, then
4.2 microseconds / step 23 / watermark 3448. All 12 injections had been
processed; the current computational-split counter was zero. With 48 injected
parcels and computational splitting disabled in the input, this observes new
IDs from the KHRT path. It does not establish final surviving population or
fuel conservation. Scripts/logs: `build/spray-a-progress.gdb`,
`build/spray-a-progress.log`, `build/spray-a-progress-2.log`. The same process
was still live at elapsed 13 min 56 s; no restart was performed.

The frozen Jet run was separately observed at 0.5556 ms / step 3169 /
next event 2778 / particle-ID watermark 11108
(`build/jet-energy-progress.log`). It remains in progress, not a completed
1.5 ms qualification.

For subsequent runs, a separate Release configuration was created at
`build/dodecane-simple-mpi-release` with the same compiler/MPI/bundle options
as the Debug build, changing only `CMAKE_BUILD_TYPE=Release`. Configure passes
(`build/dodecane-release-configure.log`); the selected spray target build is
running (`build/dodecane-release-build.log`). No Release test is counted yet.
The currently running Debug and frozen Jet binaries are unchanged.
