# Selected AMR checkpoint build isolation — 2026-09-17

Source parent: `d53c5ea` on `main`, plus this change.

The shared 3D AMR checkpoint source acquired direct fixed-H2/O2 imports when
fixed chemistry context was added. Those imports broke the selected runtime's
independent module search path. The selected MPI target also lacked the newer
`mpi_amr_checkpoint_3d_mod` source required by its application driver.

The selected AMR target now privately defines `PELEF_SELECTED_ONLY_RUNTIME`.
That build omits fixed database imports and rejects fixed chemistry context;
selected reaction/provenance context remains caller-supplied. The default core
build retains the existing fixed context implementation. The MPI source list
now includes the distributed checkpoint module. No schema, tolerance, generated
mechanism or golden hash changed.

## Reproduction

New build directory, GNU Fortran 13.3 Debug, OpenMPI 4.1.6, MPI/spray enabled:

```sh
cmake -S . -B build/selected-amr-checkpoint-isolation-debug -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran \
  -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec \
  -DMPI_mpi_LIBRARY=/usr/lib/x86_64-linux-gnu/libmpi.so \
  -DMPI_hwloc_LIBRARY=/usr/lib/x86_64-linux-gnu/libhwloc.so \
  -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_SPRAY=ON
cmake --build build/selected-amr-checkpoint-isolation-debug \
  --target test_selected_full_amr_reactive_3d test_selected_full_mpi_amr_reactive_3d \
  test_amr_reactive_3d_checkpoint -j 4
ctest --test-dir build/selected-amr-checkpoint-isolation-debug \
  -R '^(unit_amr_reactive_3d_checkpoint|regression_selected_amr_reactive_3d_(restart_(reference|run|exact)|checkpoint_stop)|regression_selected_mpi_amr_reactive_3d_(checkpoint_np2|restart_np1))$' \
  --output-on-failure
```

Initial configure/link attempts discovered Intel MPI and hwloc through the
shell environment. Explicit GNU/OpenMPI wrapper, launcher and library paths
above corrected that mismatch. The intermediate clean build exposed the MPI
checkpoint source-list omission, which is included in this fix.

`nm` confirms no fixed database module symbols in either selected AMR runtime
archive. The generated full-H2/O2 fixture executable intentionally contains its
selected generated module with the same module name; absence of that name in
the executable is therefore not an appropriate isolation criterion.

Logs: `build/selected-amr-checkpoint-isolation-debug/evidence/`.
This is focused build/restart evidence, not a full configuration matrix or
physical validation.

## Results

Final target builds passed. Focused CTest: **7/7 passed, no failures/skips,
98.21 s**. The fixed checkpoint unit covers fixed reacting context. Selected
serial restart retains the existing checkpoint and output golden hashes.
MPI saved at two ranks and restarted at one rank; direct byte comparisons of
both final CSVs against the serial continuous reference passed:

- Coarse SHA256: `a03884837909bab719f0783d23ea9150a81f2afc628fee1a01c645c1b47e0ff3`.
- Fine SHA256: `8d5a992d07b89832624e37af6ddacd510df63633707678693b9053057d4efcbe`.

Mechanism provenance preflight and project contract passed; Python AST checks
of the provenance, contract and comparison scripts passed (3/3).
`git diff --check` passed. The full target set was not rerun.
