# Serial embedded-wall transport — 2026-09-27

## Implementation

`3da8e1608da2e94b96aebd62e3bf363cac6bcc0a`, based on `5f7205bb`, connects the
existing quadratic embedded-wall heat/viscous kernel to serial EB tree SDC.
Two donor layers supply each wall stencil, including interpolated ancestor
states during recursive stepping. Composite integration excludes covered
coarse cells; heat and impulse use the same stage weights as transport.
Outputs remain private until the whole operation succeeds.

The fixed/selected serial application now accepts `embedded_wall_temperature`
and `embedded_noslip`, builds the two-layer geometry cache, caps the interval
with the wall row-sum frequency and checks wall settings on restart. This cap
also applies when automatic gas timestep selection is off, as in the MPI app.
Stationary adiabatic viscous walls exchange momentum with zero mechanical
energy transfer. No new wall model, clipping or tolerance changes are introduced.

## Completed checks

GNU Fortran 13.3.0 / OpenMPI 4.1.6, Debug, MPI/spray enabled and tests/SUNDIALS
disabled in `build/heptane-simple-mpi-debug`. The fixed core, MPI support and
fixed serial/MPI apps were rebuilt. The selected serial app compiled; its
runtime was not separately qualified here. Tested numerical source matches the
implementation commit; three comment lines were relocated after the builds.

One two-rank run of the extended `test_mpi_eb_wall_heat_3d` passes:

- Serial wall rates and frequency exactly equal distributed evaluation in both
  split ownership and all-on-rank-zero layouts. Composite heat/traction agree
  with the existing whole-fine-grid reference. Negative wall temperature
  rejects with exact rate, ledger and frequency rollback.
- Eight evolved cases (heat or adiabatic viscous traction, CTU or MOL, common
  or recursive stepping) exactly match serial/MPI gas state and temperature
  on owned patches. Ledgers agree within the unchanged `1e-12` gate.
- Cooling from 500 K toward a 400 K wall over `1e-7 s` transfers approximately
  `-4.2739859565e-6 J`; maximum normalized budget residual is `1.793e-15`.
  Tangential 2 m/s flow transfers `-5.9822404766e-11 kg m/s`, with zero wall
  energy and budget residual `3.113e-17`. The existing two wall-free MOL
  conservation cases also pass.

Eight application launches finish as expected: continuous, checkpoint stop,
serial resume, two-rank MPI resume, two wall-context rejection cases, and serial
and MPI interval-cap cases. The reacting two-level MOL fixture with both wall
options reaches `2e-10 s`; serial and MPI resumed gas CSVs exactly match
continuous serial execution. Its composite error is `2.700e-16`. Changed wall
temperature and removed no-slip each return exit 1 with the expected diagnostic
and no result CSV.

With the requested interval enlarged to `1e-6 s`, `transport_cfl=1e-14` and
checkpoint stop after one step, both executables select
`2.743960716823396e-13 s` and produce identical gas CSVs. These runs verify the
wall cap; they are stopped checkpoints, not completed `1e-6 s` trajectories.

Provenance preflight, project contract, Python syntax for the artifact driver
and `git diff --check` pass. Existing compiler/executable-stack warnings and
MPI display-authorization diagnostics remain in logs. No full CTest matrix,
installed-binary audit, new external parity or physical validation is claimed.

## Artifacts

`build/serial-eb-wall/commands.txt` records commands. `build.log`,
`selected-build.log`, `check.log` and `app.log` retain outcomes. `run_app.py`
uses the preceding `build/serial-mol/app/continuous/tree.nml` fixture;
`app/result.json`, `app/runs.json` and per-run inputs preserve commands,
input hashes, exit codes, durations and binary hashes. Repetition needs a fresh
application output directory.

| Executable | SHA-256 |
| --- | --- |
| Fixed serial | `b8c2059fc67d310fb5158653911556623b3f881c3b0320569541a5bbaa0b9d3c` |
| Fixed MPI | `f5116d6618f362e81349240f34d7ff602e965a443b879b31e12f444b2aa4f023` |
| Selected serial, build only | `15e16775d93f337a08fb3108d096a54ea55c2443f8caa3373d145186d9d24c7a` |
| Extended wall test | `42f61afc13c4a01fdf897fcf741cb34ec3e8ab50e92c9127e03292f2ac277634` |
