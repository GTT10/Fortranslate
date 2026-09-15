# Mixed-film trace depletion — 2026-09-16

Tested parent: `37881f7` plus this change. Source/test hashes and build/test/preflight
logs are retained in `build/mixture-film-trace-checks/`.

Actively evaporating film components at or below 1e-12 of current film mass, or
whose predicted residual reaches that threshold, now transfer the full residual
through the conservative exchange. Frozen and grid updates enable this policy;
the supplied-exchange API defaults to exact supplied amounts unless its new optional
flag is enabled. Suppressed components and zero intervals are unchanged. Relative
depletion limits exclude already trace components, preventing their shrinking
inventories from imposing progressively smaller steps. Bulk overdraw still rejects.

CLI mixed-wall contexts add `mixed-film-v2-trace`; the prior mixed-wall evolution
therefore cannot resume silently. API callers must include the public evolution ID.
File layout and pure-wall/particle-only mixture context are unchanged.

GNU Fortran 13.3.0 Debug, MPI disabled, bounds/FPE checks enabled:

```sh
cmake --build build/port-spray-debug --target test_spray_core pelef_spray_3d --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(mixture_wall_rates|mixture_wall_exchange)$' --output-on-failure
python3 tests/spray/test_spray_application.py --exe build/port-spray-debug/pelef_spray_3d SprayApplication.test_mixed_wall_cli_restart
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Two CTest entries passed in 0.03 s total and one selected CLI method passed in
4.477 s, no failures/skips. New cases verify complete active trace transfer to its
gas species, ambient-suppressed trace retention and zero-interval preservation.
Existing film exchange/conservation and mixed-wall CLI restart cases pass.
Provenance, contract and whitespace checks passed; no Python source changed.
No full suite, configuration matrix or physical validation was run.
