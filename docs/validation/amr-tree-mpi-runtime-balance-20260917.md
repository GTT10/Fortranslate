# MPI runtime work sampling and rebalancing — 2026-09-17

## Source and scope

Base `af07cf9822d833c1afdb0b0bcdc0bb2bc27b102e`, local `main`.
Adds opt-in per-patch kernel timing and step-cadence rebalancing to fixed/selected
MPI tree applications. Numerical formulas, mechanisms, tolerances and golden
hashes are unchanged. Default cadence zero leaves sampling disabled.

SHA-256:

- `src/parallel/mpi_amr_tree_work_mod.F90`:
  `adc91a928a958b403e99fd7fcaf5a2b27ec0390b625da6611b706b481f1b106b`
- `src/driver/mpi_amr_patch_tree_application_3d_mod.F90`:
  `34f139c4991a36bfe57baffcbc4d83fe7348d332d6572f04d9b783e0512e76f6`
- `src/parallel/mpi_amr_patch_tree_chemistry_3d_mod.F90`:
  `5bc791e893ad8026178c7a5c6ec90ec7a90754c327aa7d3e8cac659450d0c6e1`
- `src/parallel/mpi_amr_patch_tree_hydro_3d_mod.F90`:
  `375e13286dbd1d9afc106df2a2a59967fd7e0b9115f4fec4d7e7283637834bdb`
- `src/parallel/mpi_amr_patch_tree_transport_3d_mod.F90`:
  `0d83d548deb4cf18c95d92a74fb2eec8b7fcfce84f46d6f06c113e000c8bcf54`
- `src/parallel/mpi_amr_patch_tree_hydro_level_3d_mod.F90`:
  `7661441e42ce663c9703b058e8ee23b3fe6ea2f5217b88d14d3f4590cafc3ad8`
- `src/parallel/mpi_amr_patch_tree_transport_level_3d_mod.F90`:
  `e7c624e0206c7471c19664aee481d6417e2ceaaa169b8404d87daee41cb263b6`
- `tools/check_mpi_amr_tree_restart_3d.py`:
  `ffdc543e8634eda7b85f782f1ea67274bf94d8b9478e494db393a3d0556b03b4`

## Commands

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6; existing matching MPI build.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target pelef_mpi_amr_tree_reactive_3d test_selected_elementary_mpi_amr_tree_3d -j 4
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target pelef_mpi_amr_tree_reactive_3d test_selected_elementary_mpi_amr_tree_3d pelef_amr_tree_reactive_3d test_selected_elementary_amr_tree_3d -j 4
python3 -m py_compile tools/check_mpi_amr_tree_restart_3d.py
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^regression_mpi_amr_tree_(recursive_)?runtime_balance_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

## Completed results

Builds passed. CTest **2/2 passed**, no failures/skips, 108.06 s total:
common 72.35 s and recursive 35.68 s. Provenance, Python syntax, project contract
and whitespace checks passed. Existing compiler warnings remain. Runtime mechanism
coverage uses fixed/selected elementary fixtures, not a broad detailed-mechanism
matrix. Common continuous runs complete five steps and migrate at step one for
both mechanisms. Recursive runs complete two steps and retain their initial owners.

## Coverage and limits

The checker enables every-step sampling with zero minimum gain on MPI runs.
It compares full-split fixed/selected common and recursive results to serial bytes,
including budget CSV and checkpoints. Two-rank saves continue on three ranks with
rebalancing disabled, on one rank with it enabled, and in serial. Negative cadence
and gain 1 reject before output changes. Existing checkpoint/context guards remain.
Actual migration is observed in common-step fixed/selected continuous logs; tests
do not require a particular noisy timing outcome or claim a measured speedup.

Timers include owner-local flux/update/packing and chemistry solves. Ghost
communication, restriction, standalone EOS recovery, timestep selection and
migration overhead are excluded. Failed-attempt samples are discarded by code;
this regression does not force a retry. The process-local sampling scope is not
thread-safe or reentrant. Dynamic geometry combined with runtime rebalancing,
full configuration matrices, installed audits, production scaling, external solver
parity and physical validation were not exercised.

Logs: `build/amr-tree-mpi-runtime-balance-checks/{build,rebuild,ctest}.log`.
Artifacts: `build/selected-tree-mpi-debug/tests/mpi_amr_tree_runtime_balance_3d/`
and `mpi_amr_tree_recursive_runtime_balance_3d/` in the same parent directory.
