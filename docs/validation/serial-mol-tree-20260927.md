# Serial EB tree MOL connection — 2026-09-27

## Implementation

Source `93723bf31daddf0243ae534ea89b49d96de4563c`, based on `3b9dedb1`,
connects the existing MOL reconstruction, interface flux and EB centroid
operators to serial tree common and recursive evolution. A three-cell halo
adapter evaluates each patch, including ancestor endpoint interpolation.
Hydro and diffusion rates/fluxes are combined at the old and new stages;
nonreacting steps use one trapezoidal correction and reacting iterations
re-evaluate the preceding reacted endpoint. Embedded slip-wall pressure
impulses use the same stage weights. Existing redistribution, reflux and
transactional publication remain in use.

Serial fixed/selected applications now accept `mol_mode` and
`mol_lax_friedrichs` with `time_integrator='sdc'`. Existing MOL/LF checkpoint
identity applies to serial saves and continuation. CTU remains the default.
This does not add serial embedded-wall heat/no-slip transport or serial spray,
and does not repair the known nested LF negative-fuel case. Reconstruction
retains the previously documented positive-face extension to the pinned source.

## Completed focused checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Debug, existing
`build/heptane-simple-mpi-debug`; MPI/spray enabled, tests/SUNDIALS disabled.
The fixed core, MPI support and fixed serial/MPI applications were rebuilt.
The selected serial application also compiled; selected runtime behavior was
not independently exercised in this increment.

One two-rank invocation of `test_mpi_eb_tree_full_3d mol` passes all four
variants below. The fixture has five patches, three levels, a planar cut at
z=0.31, nonuniform gas, forcing and a `1e-10 s` root interval. Gas state and
temperature match serial evaluation **exactly** on every owned patch.
Ledger differences remain below the existing `2e-12` gate, and normalized
composite closure remains below the existing `3e-14` gate.

| Variant | Stepping | Flux | Chemistry | Diffusion | Redistribution | Composite residual |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | Common | Acoustic | Off | Off | None | `4.1725113643490740e-15` |
| 2 | Common | LF | On | All | None | `4.3482811143688856e-15` |
| 3 | Recursive | LF | Off | All | Flux | `4.1011907183895860e-15` |
| 4 | Recursive | Acoustic | On | All | State | `4.3410205888251901e-15` |

One invocation of `test_serial_sdc_species_3d` passes four source cases:
negative-species rejection with complete rollback and positive trace-species
acceptance, each with CTU and MOL. The registered new MPI CTest is
`unit_mpi_serial_mol_tree_3d_np2`; these were direct executable runs, not a
CTest matrix.

Six application launches complete as expected: serial two-step continuous,
checkpoint stop, serial resume, two-rank MPI resume, and expected rejection of
CTU/LF mode changes. The two-level reacting MOL fixture uses a 2³ root,
ratio-two child and a `2e-10 s` final time. Both resumed gas CSVs are byte-for-byte
identical to continuous serial execution. The two deliberate rejections return
exit 1 and produce no result CSV. These are internal implementation/restart
checks, not external PeleC parity or physical validation.

Provenance preflight, project contract, artifact-driver Python syntax and
`git diff --check` pass. An initial development compile forwarded a MOL keyword
to the reactor-only API; that call was corrected before the successful build.
An initial manually linked MPI test inherited a foreign library search path
and could not load `libmpi.so.12`; rebuilding with `LD_LIBRARY_PATH` and
`LIBRARY_PATH` unset produced the successful OpenMPI run. Existing compiler,
executable-stack and display-authorization diagnostics remain in logs.

## Artifacts

`build/serial-mol/commands.txt` contains commands; `build.log`,
`selected-build.log`, `parity.log` and `species.log` contain outputs.
`run_app.py` and `app/result.json` retain the six application inputs/commands,
exit codes, durations and binary hashes. `run_app.py` requires a new output
directory when repeated.

| Executable | SHA-256 |
| --- | --- |
| Fixed serial app | `133be84a693426923f2c04f124cea6458d71caff18fee08bb9767a1658d0081b` |
| Fixed MPI app | `7941873db0c197494c783d74b6ec72202d4d7c14e5ae571abb56c40f8c150ec6` |
| Selected serial app (build only) | `cb02346386355a19536eddc9ff691da6a34d4b138b9ecd653fb1ff4bad837b95` |
| MPI parity test | `d86d2b556f090c461aa2a5dba8a0ffbc6bf1a306875548d7f51fda6f4f416f6c` |
| Species test | `7bf24a4cb37c1f688ac3c0719874b7c57a3101c06ad3a62abe689f5e4a703cc4` |
