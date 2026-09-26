# MPI EB tree SDC checkpoints — 2026-09-27

## Source and format

Implementation `a73a1c10`, based on `80b5168b` on local `main`.
PeleC pin remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

MPI writers select V4 for SDC and retain V2 for split stepping. V4 adds the
SDC iteration count, chemistry/source presence, external-force configuration
and cumulative forcing inventory after the existing tree metadata. Indexed
patch payloads contain a legacy gas record followed, when chemistry is enabled,
by `SDC_REACTION_SOURCE` and the signed reaction-estimate rows. Each owner writes
its own payload; source arrays are not gathered to rank 0. The outer envelope
binds SDC controls while embedded gas records retain mechanism/configuration
checks through a private split configuration.

MPI readers use direct V4 offsets, retain V1/V2 support and now also read serial
V3 SDC records. V3 source offsets are found by a root scan, then each owner reads
its own source rows. Serial readers accept V4 as well as prior formats. Forcing
participates in composite inventory closure. Read state, temperature and the
checkpoint record publish only after all payload, envelope, finite-source,
inventory and restriction-consistency checks succeed collectively.

## Completed focused checks

GNU Fortran 13.3.0, configured `/usr/bin/mpifort` and `/usr/bin/mpiexec`
(OpenMPI/OpenRTE 4.1.6), Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`. Existing compiler warnings remain.

- Mechanism provenance preflight passed for `mechanisms/h2o2_full.json`.
- The first build failed on one Fortran line exceeding 132 characters. It was
  wrapped; no test ran against that failed build. The subsequent build passed.
- CTest: 2/2 passed, 1.43 s total: existing `unit_eb_tree_checkpoint_3d` (0.51 s)
  and new `unit_mpi_eb_tree_sdc_checkpoint_3d_2` (0.87 s).
- New synthetic checkpoint fixture uses three levels, five patches and two
  writer ranks. Both chemistry-enabled signed reaction estimates and
  chemistry-disabled absence are exercised.
- Round-robin and all-on-rank-0 reads preserve gas, temperature, reaction
  estimates and forcing exactly, with non-owned arrays unallocated.
- Each process separately reads the two-rank file with `MPI_COMM_SELF`, proving
  the one-rank reader path without changing the file. Gas and reaction estimates
  match exactly. This is checkpoint data interchange, not a continued trajectory.
- Serial readers accept MPI V4, and MPI readers accept serial V3; tested source
  and gas values are exact.
- Mismatched external-force configuration rejects without publishing. A NaN row
  in the last patch's reaction payload also rejects and preserves prior gas,
  temperature, source estimates and record time.
- Python syntax, project contract and `git diff --check` passed.

Commands from repository root:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_sdc_checkpoint_3d test_eb_tree_checkpoint_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_(mpi_eb_tree_sdc_checkpoint_3d_2|eb_tree_checkpoint_3d)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Logs: `build/sdc-mpi-eb-checkpoint-checks/`. Test checkpoint files are written in
`build/spray-amr-reactive-debug/tests/`.
No MPI EB SDC application continuation, selected build, full configuration
matrix, installed-binary audit, regrid, recursive stepping or physical
validation was run. Application policy-header routing remains separate work.
