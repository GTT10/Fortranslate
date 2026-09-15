# Computational parcel splitting — 2026-09-15

Source: parent `43abdba` plus this commit's parcel/breakup/coupling/checkpoint
and application code, checker, tests and docs. GNU Fortran 13.3 Debug,
`build/port-spray-debug`.

## Implemented boundary

ETAB can cap represented droplet multiplicity by partitioning a parcel into
multiple computational parcels. It preserves droplet size, velocity,
temperature, position and deformation, and divides represented weight with a
remainder correction. No transverse velocity or stochastic dispersion is added.
This is computational population management, not a new physical breakup law.

All child counts, storage capacity and integer ID bounds are checked before
replacing arrays. Existing IDs remain stable; children append in parent order
and store their immediate parent ID. The coupled step keeps candidates until
splitting succeeds, so failure preserves gas, temperature, parcels and ledger.
Scheduled injections use IDs above the current maximum. Restart counts original
injection roots, not all descendants, against the injection schedule.

V5/V6 checkpoints persist parent IDs as strict integer records in addition to
V3/V4's deformation and optional physical-boundary ledger. Reads require each
parent to exist at an earlier ID, using binary search over sorted predecessors.
Writers reject silently dropping lineage into an older format. Limits are bound
into the splitting context. Existing non-splitting formats remain supported.

## Focused evidence

- `spray_split`: four parcels of multiplicity 10 become sixteen parcels with
  weights <=3. IDs are canonical, roots and children are distinguished, and
  aggregate mass/momentum/energy are conserved. An insufficient storage limit
  rejects without changing allocation, IDs or weights.
- Initial `spray_split_restart`, `spray_split`, existing
  `spray_integration_checkpoint` and `spray_application_cli`: 4/4 passed in
  18.06 seconds.
- The restart checker was extended with two scheduled injections following the
  checkpoint. Periodic and physical-wall runs preserve gas/parcel CSVs and the
  final history row exactly across stop/restart; the extended test passed in
  10.95 seconds. The final periodic run has 16 computational parcels, 12 injection
  roots and maximum multiplicity 9,441,701.201762173 below the 1e7 cap. Total
  represented liquid mass is 1.1e-10 kg within relative 1e-12.
- An additional direct replay changed the first checkpoint parcel's parent to
  itself. The reader rejected the invalid parent before producing outputs.
- Project contract, pinned mechanism provenance, Python syntax and diff
  whitespace checks pass. No full compiler/backend/MPI matrix was run.

Scheduled-injection artifacts and the invalid-parent probe are in
`build/port-spray-debug/tests/spray_split_restart/etab-mblb3b20/`.

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d test_spray_core test_spray_integration --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(split|split_restart|integration_checkpoint|application_cli)$' --output-on-failure
python3 tools/check_project_contract.py --source .
python3 tools/check_mechanism_provenance.py --bundle build/port-methanol-reference/ffcm1_red.json --source build/port-methanol-reference/ffcm1_red.yaml
python3 -m py_compile tools/check_spray_etab.py
```

This fixture disables drag/evaporation and isolates breakup, splitting,
scheduled injection and persistence. Partition-count independence of nonlinear
sequential gas feedback, transverse child dispersion, particle MPI/AMR/EB and
experimental breakup validation remain open. Final allocation guards were
compiled after the runtime tests.
