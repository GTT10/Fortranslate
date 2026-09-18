# EB impact child ID ordering — 2026-09-18

Test implementation `1fe8fdbbbf5da544e9b0baa68f04fc81e7496e9c`, based on
`c410c20`; solver behavior comes from `026dac6`.

Inspection confirmed `impact_mpi_level_set` delegates child construction to
the shared `resolve_mpi_wall_impacts` routine. No second implementation was
needed. The EB fixture now places parent IDs in reverse owner-rank order and
checks exact parent-ordered child IDs and lineage after oblique level-set
contact and wet splash.

## Completed verification

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug, existing fixed configuration.
`test_spray_mpi_eb_impact` rebuilt. CTest **2/2 passed**, zero failures/skips:
`spray_mpi_eb_impact_2` (0.46 s) and `spray_mpi_eb_impact_4` (0.47 s), 0.97 s total.

Every original parent produces four children with exact ID formula
`nranks + 4*(parent_id-1) + sibling`, and the global watermark is `5*nranks`.
The existing contact-fraction, outgoing-velocity and global inventory checks
pass, with inventory tolerance `3e-12` unchanged. Existing empty-owner,
capacity, geometry disagreement and remote thick-film rejection paths also
pass. This verifies the EB wrapper's use of the shared ordering; it is not
full reactive EB restart, arbitrary topology/event ordering or physical impact
validation. No numerical algorithm changed in this increment.

Provenance, Python syntax, project contract and whitespace checks passed.
No failed builds/tests, full matrix, selected-runtime run, external solver
parity or scaling qualification.

## Commands and artifacts

Logs: `build/spray-eb-id-checks/{provenance,build,ctest}.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_impact --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_impact_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
