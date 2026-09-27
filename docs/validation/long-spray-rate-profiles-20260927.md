# Full-length external spray rate profiles — 2026-09-27

Base: `549bbbeb`. The MPI EB tree spray input now reads external rate tables
into allocated per-port profiles. The shared reader/integrator supports up to
16,384 points; fixed-size callers retain their supplied capacity. Inline
namelist profiles remain limited to 64 points. Invalid input preserves the
previous profile. No rate resampling, integration-rule or tolerance change.

The jet preparation tool now exports every one of the 7,502 pinned ROI knots
in SI units. Its two-event schedule remains unchanged; accepting the complete
table does not claim a full-duration jet run.

## Restart identity

For short external profiles, values still populate the existing namelist
arrays before context serialization. Longer tables append a tagged port/count
record and every time/rate pair at 18-digit scientific precision. Thus even a
change outside the scheduled injection window changes the restart identity.
The MPI EB particle context reader allocates from the expected context length,
with a bounded 16,008,192-character maximum; application profile context is
bounded at 16,000,000 characters. No checksum replaces exact numerical matching.
Legacy short-context layout is retained; an old-binary checkpoint migration
was not independently rerun in this increment.

## Completed checks

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6, selected decane Simple bundle,
MPI/spray ON, tests/SUNDIALS OFF:

```sh
python3 tools/check_mechanism_provenance.py --bundle build/decane-3sp-simple-reference/decane_3sp.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/decane-simple-openmpi-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
python3 tools/prepare_jet_reference.py --bundle build/decane-3sp-simple-reference/decane_3sp.json --reference-source build/pelec-pinned-hit-source --output build/decane-jet-long-rate
```

Build logs: `build/long-rate-build.log` and `build/long-rate-final-build.log`.
Both builds complete; the second includes removal of the remaining application
64-point check. Existing real-comparison and array-temporary warnings remain.

Two focused unit programs were compiled directly with `gfortran -g -fcheck=all
-ffpe-trap=invalid,zero,overflow`, using `precision_mod.F90` and
`spray_injection_profile_mod.F90`:

- `tests/spray/test_spray_rate_file.F90`: existing parsing, capacity and
  transactional rejection checks pass.
- `tests/spray/test_spray_long_rate_file.F90`: 10,001-point linear profile,
  independent analytical partial-interval integral, and nonmonotone-file
  rejection without changing the prior profile pass. This test is registered
  with CTest; the current build has tests OFF, so execution was direct.

Artifacts: `build/long-rate-unit/`.

The three-step single-rank jet run completes using the complete table. Gas and
particle CSVs are byte-identical to `build/decane-jet-connection/`, whose
four-knot file covered the same interval. The run command is the case guide's
MPI command from `build/decane-jet-long-rate/`.

A checkpoint after step one resumes on two ranks with gas CSV byte-identical
to the continuous run. Changing rate row 7,000, outside the simulation interval,
rejects restart with `EB spray context mismatch` and no final gas CSV. Exact
commands/results: `build/long-rate-restart-check/{commands.json,result.json}`;
driver: `build/check_long_rate_restart.py`.

Python syntax, project contract and diff checks pass. No full matrix, full
3.8 ms jet, event-schedule capacity expansion, or physical validation was run.
