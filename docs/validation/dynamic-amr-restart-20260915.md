# Serial dynamic 3D AMR restart — 2026-09-15

Source: parent `0bc12fe2597320799f0202471f6068ba549f85ff` plus the
checkpoint/driver/configuration/tests/docs changes in this commit. GNU Fortran
13.3, existing Debug `build/port-spray-debug` and Release
`build/port-spray-cvode`; this is focused local evidence, not a full matrix or
physical validation.

Schemas 6--9 add a versioned post-step regrid phase, exact tagging configuration,
original input bounds, active patch bounds and completed regrid count to the
existing fixed/selected and transport context. The allocatable reader validates
into candidates before committing geometry and arrays. Static schemas 2--5
remain supported. MPI dynamic regridding and fixed-mechanism chemistry restart
remain rejected.

## Completed checks

- Debug `unit_amr_regrid_3d` and `unit_amr_reactive_3d_checkpoint`: 2/2 passed.
  The checkpoint test also exercises schema 9 selected/transport roundtrip and
  rejected-read transaction preservation. Existing static schema tests pass.
- Release `tools/check_dynamic_amr_3d_restart.py`: uninterrupted advection versus
  checkpoint at coarse step 5 and restart to 2 microseconds. The patch moves at
  step 0 and expands at step 4. Both outputs report two completed regrids.
  Coarse and fine CSVs match byte for byte. Changed tagging threshold and a
  truncated checkpoint are rejected.
- Debug `amr_dynamic_restart_3d`: 1/1 passed in 67.44 seconds after the
  candidate-allocation and transaction changes; artifacts are in
  `build/port-spray-debug/tests/amr_dynamic_restart_3d/`.
- Project contract, Python syntax, source provenance and `git diff --check` pass.
  Pinned FFCM1_Red YAML SHA256:
  `ddbba8c97d929530cfb213da45485cca4f979dccf97c1063721bd855d7b3f259`.

Release artifacts: `build/port-spray-cvode/dynamic-restart/` (inputs, logs,
checkpoint, reference/resumed CSVs). SHA256:

- Coarse: `b5cd4361685a1d4cea13331d095d1a7f92bdd494ccb72695384c8907732b8808`
- Fine: `aa36eb4ea6ee2bf980b0defe42b75ddff2242155834b9b0e26d3af924d677687`

Commands:

```sh
cmake --build build/port-spray-debug --target pelef_amr_reactive_3d test_amr_regrid_3d test_amr_reactive_3d_checkpoint --parallel 8
ctest --test-dir build/port-spray-debug -R '^unit_amr_(regrid_3d|reactive_3d_checkpoint)$' --output-on-failure
python3 tools/check_dynamic_amr_3d_restart.py --executable build/port-spray-cvode/pelef_amr_reactive_3d --work-dir build/port-spray-cvode/dynamic-restart
ctest --test-dir build/port-spray-debug -R '^amr_dynamic_restart_3d$' --output-on-failure
python3 -m py_compile tools/check_dynamic_amr_3d_restart.py
python3 tools/check_project_contract.py --source .
python3 tools/check_mechanism_provenance.py --bundle build/port-methanol-reference/ffcm1_red.json --source build/port-methanol-reference/ffcm1_red.yaml
```

The initial CTest filter omitted the `unit_` prefix and matched zero tests;
the corrected filter above ran both tests. An initial project-contract command
omitted `--source` and printed usage; the corrected command passed. No full
compiler/backend/MPI matrix or mechanism regeneration was run in this increment.
