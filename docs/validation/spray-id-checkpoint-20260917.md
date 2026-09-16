# Persistent spray ID high-water mark — 2026-09-17

Parent `7e30bcf` plus this change. Build/test/preflight logs and source hashes:
`build/spray-id-checkpoint-checks/`.

Optional `last_id` on checkpoint read/write selects `_ID_V1` and stores the
64-bit reservation watermark as a strict integer record. It must cover all live
IDs and can exceed them, retaining reservations for removed particles. Readers
commit the watermark only after complete parsing and state validation. Omission
retains the existing checkpoint format. MPI collection/scatter and CLI restart
are not yet connected by this API increment.

GNU Fortran 13.3.0 Debug, MPI/SUNDIALS disabled in existing `build/port-spray-debug`:

```sh
cmake --build build/port-spray-debug --target test_spray_integration --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_integration_checkpoint$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

1/1 focused checkpoint entry passes (0.06 s total), no failures/skips. Existing
legacy and KH timer roundtrips pass. New checks preserve `huge(int64)-7` exactly,
require explicit ID mode, reject a watermark below live IDs, and reject an
out-of-range integer token without changing gas/temperature or the caller's
watermark. Provenance, project contract and whitespace checks pass; no Python
source changed. No full matrix or changed-rank MPI restart is claimed.
