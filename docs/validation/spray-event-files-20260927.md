# External spray event schedules — 2026-09-27

Base: `defa8ad4`. The MPI EB tree spray application now accepts an external
schedule with up to 65,536 events, removing the 1,024-event inline limit for
that input route. The four columns are begin time, end time, integer port and
scale. Parsing is transactional and requires finite bounded values,
nonnegative ordered intervals, nondecreasing starts and ports 1–16. The
application additionally validates the configured port count and profile
coverage. Inline input remains available and cannot conflict with a file.

Allocated schedules supply timestep clipping, interval mass/velocity
preparation, event dispatch and restart event checks. Exact numerical rows
are appended to checkpoint context under `EVENT_FILE_V1`; filenames do not
replace numerical identity. The existing combined context capacity still
applies. Inline contexts retain their existing representation.

The jet preparation tool writes a schedule file and accepts `--events`.
All 7,499 active intervals through 1.5 ms were prepared in
`build/jet-events-full-prepared/`; **that long run was not launched**. The
default remains two active intervals. This does not establish full Spray-Jet
trajectory parity or remove the documented mesh/material adaptations.

## Completed checks

GNU Fortran 13.3.0 Debug, OpenMPI 4.1.6, selected decane Simple bundle,
MPI/spray ON, tests/SUNDIALS OFF:

```sh
python3 tools/check_mechanism_provenance.py --bundle build/decane-3sp-simple-reference/decane_3sp.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/decane-simple-openmpi-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
python3 build/check_spray_event_application.py
python3 -m py_compile tools/prepare_jet_reference.py build/check_spray_event_application.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The build passes (`build/spray-event-file-build.log`). The new
`tests/spray/test_spray_event_file.F90` is registered with CTest and was
compiled/executed directly with `gfortran -g -fcheck=all
-ffpe-trap=invalid,zero,overflow`, linking the precision/profile modules.
It reads 2,048 events and rejects reversed intervals, invalid/noninteger
ports, excess fields and unsorted starts while retaining the previous table.
Artifacts: `build/spray-event-unit/`. No CTest matrix was run.

The application fixture processes **1,100 events** in three gas steps. Two
events have positive mass (including event 1,100); 1,098 have zero scale.
This deliberately tests indexing beyond the old bound without a long gas
simulation. Gas and particle CSVs are byte-identical to the preceding inline
two-event baseline. It is not evidence for 1,100 positive-mass injections.

A step-one checkpoint resumes on two ranks with identical final gas output.
Changing the scale of event 1,100 rejects restart with `EB spray context
mismatch`, without producing the final gas CSV. Exact commands and results:
`build/spray-event-application/{commands.json,result.json}`.

The preceding implementation's inline-schedule checkpoint also resumes on
two ranks with gas output byte-identical to the original continuous run
(`build/spray-events-legacy-restart/run.log`). Its first launch used the wrong
working directory and failed to locate the executable; it is retained in
`initial-launch-failure.log` and is not counted as a solver run.

Provenance, Python syntax, contract and diff checks pass. No full-duration jet,
full configuration matrix, or new physical validation was performed.
