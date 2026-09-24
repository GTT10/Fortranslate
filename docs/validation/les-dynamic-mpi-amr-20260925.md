# Periodic Dynamic MPI AMR operator wiring — 2026-09-25

Implementation `e457cb9f7ea2ebebfbfd19c302c4351a05166ff4`, based on `3bdb6bc`.

The MPI AMR LES timestep selector and SSPRK2 operator now accept periodic
`model='dynamic'` and optional test filter/model ratio, grid ratio and harmonic
face interpolation, with the serial API defaults. Effective filter strings and
numeric options join rank agreement before donor communication. The agreed
filter determines halo width. Each owner evaluates dynamic faces and a padded
diffusion bound, then uses the existing distributed restriction/reflux and
transactional state/temperature/diagnostic publication. Static-only options
remain compatible; dynamic controls on static modes reject.

## Evidence boundary

GNU Fortran 13.3.0 / Open MPI 4.1.6, Ninja Debug, existing configurations.
Fixed `test_les_mpi_amr_tree` and selected H2O2 `pelef_selected_spray_tree_runtime`
builds passed. Fixed test was rebuilt after adding the filter disagreement check.
Provenance preflight, Python syntax, project contract and whitespace passed.

**No MPI runtime test was executed for this change.** The earlier launcher and
CTest PIDs 2018713/2018710 were still present at this turn's initial observation.
No repeated launch was attempted. New test cases are compiled, not passed:
serial/MPI dynamic timestep/state/temperature/ledger comparison for root-only and
nested trees, empty-rank and distributed ownership, oversized-interval rollback,
and rank-disagreeing dynamic filter ratio rejection. Existing tests already
cover generic model/CI/CFL agreement and rollback.

This is implementation/build evidence only. Reacting-spray/application dispatch,
Dynamic nonperiodic boundaries and LES subcycling remain separate work. No
physical validation, external solver comparison, scaling or full matrix ran.

## Commands and artifacts

Logs: `build/les-dynamic-mpi-amr-checks/` (`provenance.log`, `build.log`,
`build-final.log`, `build-selected.log`, `contract.log`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_mpi_amr_tree --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Recheck the existing launcher before starting the focused 2-rank MPI AMR LES
test. A later runtime result must be recorded separately from this build-only
record.
