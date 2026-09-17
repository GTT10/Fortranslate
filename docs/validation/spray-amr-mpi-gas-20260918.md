# Distributed AMR spray gas feedback — 2026-09-18

Implementation: `41e3f4d02d586042b7045b4863c160884fe4837f`, based on
`28dfc84`. `deposit_mpi_amr_tree_spray_gas` validates the initial owner-local
gas, deposits into private conserved fields, then recovers EOS temperatures
after covered-parent averaging. Conserved fields and temperatures publish
together only after collective success. Empty feedback validates and preserves
both original fields. The caller supplies identical ordered species tables on
all ranks; request and receive capacities remain available.

## Completed focused verification

Fresh GNU Fortran 13.3.0 Debug/Ninja builds used `/usr/bin/gfortran` and matching
Open MPI `/usr/bin/mpifort`, `/usr/bin/mpiexec`. Fixed tests were enabled;
the selected H2/O2 library build disabled tests. No install prefix was used.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/spray-amr-mpi-gas-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_SPRAY=ON -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=ON
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-mpi-gas-debug --target test_spray_mpi_amr_gas --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-mpi-gas-debug -R '^spray_mpi_amr_gas_(2|4)$' --output-on-failure --output-log build/spray-amr-mpi-gas-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/spray-amr-mpi-gas-selected-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_SPRAY=ON -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=OFF -DPELEF_MECHANISM_BUNDLE=mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-mpi-gas-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, no failures/skips: two ranks **0.53 s**, four ranks **0.55 s**,
**1.11 s** total. The three-level fixture compares all owned conserved cells and
temperatures exactly with serial feedback in matching source order, including
averaged parents. Its serial reference checks composite budgets within 3e-13
of the per-component inventory/contribution scale. Tests cover changed owners,
globally empty feedback, and collective rollback for negative density, impossible
energy, negative species and receive-capacity rejection. Invalid initial gas
and missing owner temperature storage also reject.

The selected library compiled; its Ninja dependency edge contains selected
runtimes without `pelef_core`. This is compile/dependency evidence, not selected
executable qualification. Provenance, Python syntax, project contract and
whitespace checks passed. No build/test attempt failed; compiler warnings remain.
No Python source changed. Allocation-failure injection, MPI fault recovery,
coupled particle stepping, scaling, AMR convergence and physical validation
were not performed.

Logs: `build/spray-amr-mpi-gas-checks/{configure,build,ctest,configure-selected,build-selected}.log`.
Executable: `build/spray-amr-mpi-gas-debug/tests/test_spray_mpi_amr_gas`.
