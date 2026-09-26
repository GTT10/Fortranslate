# Recursive spray regrid and changed-rank continuation — 2026-09-27

Regression `c313a7274a4094bf06b741b2982719f124bf713c`, based on `6f1de969`,
local main, version 0.246.0. No numerical solver changes were required.

## Verified path

The fixed MPI EB spray application advances two 10 ns recursive SDC intervals
with chemistry, molecular transport and gas forcing. At 10 ns, the refined
half-domain changes from left to right and the physical particles and SDC
reaction fields are remapped. A two-rank run checkpoints immediately after
that regrid with the second injection event still pending. Four ranks restore
that checkpoint, inject the pending parcels once, and finish at 20 ns.

The checker compares against an uninterrupted four-rank run:

- Gas/temperature CSV rows compare exactly (216 exposed cells in this fixture).
- ID-sorted particle CSV rows compare exactly, with final IDs [1,2,3,4].
- All initial/final, wall, boundary, chemical, injected and force budget columns
  satisfy the existing 3e-12 comparison criterion.
- The stopped run reports one processed event; both completed runs report two.
- Both uninterrupted and stopped runs report the scheduled regrid.
- The gas checkpoint has two `SDC_REACTION_SOURCE` sections. Both gas and
  particle payloads contain the recursive physical-source synchronization tag.
- Attempting common-step continuation rejects with `Restart subcycling mode
  mismatch` before result output. The original three checkpoint payloads remain
  byte-identical after both continuation and rejection.

This exercises one two-level topology transition and changed-rank continuation.
It does not establish three-level coupled correctness, arbitrary rebalancing,
selected-mechanism restart, temporal order, external solver parity or physical
validation. The conservative synchronization remains a PeleF extension described
in the [source synchronization record](spray-recursive-source-sync-20260927.md).

## Completed checks

Artifacts: `build/spray-recursive-restart-regrid-checks/`.
Final case directory:
`build/spray-amr-reactive-debug/spray-sdc-recursive-restart-regrid/recursive_restart_9o6z4m72/`.
Existing fixed-mechanism Ninja Debug build, GNU Fortran 13.3.0,
`/usr/bin/mpifort` and `/usr/bin/mpiexec` OpenMPI 4.1.6.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc_recursive_restart_regrid$' --output-on-failure
python3 -m py_compile tests/application/check_recursive_eb_spray_restart.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **1/1 passed**, 71.09 s (71.14 s total), no skips. Four MPI launches
comprise the single test: uninterrupted, stopped, continued and expected mode
rejection. Provenance, build, Python syntax, project contract and whitespace
checks passed. The solver executable was unchanged from the preceding app work.

Two preceding checker failures are retained:

1. 52.32 s: the checker searched the manifest for the mode tag. The actual tag
   resides in the gas and particle payloads; the assertion now checks both.
2. 70.90 s: the checker expected a later spray-context error, but the existing
   subcycling-mode check correctly rejected earlier. The assertion now names
   that actual rejection. Independent inspection of these generated outputs
   already showed exact gas, particle and budget CSV equality.

No physical run failed and no numerical/comparison tolerance changed.
`ctest.log`, `LastTest-first.log`, `ctest-message-check.log` and
`LastTest-message-check.log` retain those checker failures. `ctest-final.log`
and `LastTest-final.log` contain the final pass. No full configuration matrix
or unrelated regression suite was run.
