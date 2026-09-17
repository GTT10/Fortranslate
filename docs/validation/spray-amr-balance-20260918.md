# Gas-work-based MPI AMR spray balancing — 2026-09-18

Implementation: `fcad48e6f85398622d3adf85b12234d17357b13c`, based on
`ffe710b`. `rebalance_mpi_amr_spray` wraps the existing measured-cost gas
balancer with private owner/U/T and particle/component candidates. Ownership
and both populations commit only after final particle migration succeeds.
No-benefit or zero-cost decisions preserve the original ownership/state.

The fixed/selected application samples existing patch gas timers at configured
accepted-step intervals, restarting the sample on each timestep retry. Only
accepted-attempt costs are used. Rebalancing precedes the next adaptation,
injection and checkpoint processing. The checkpoint remains rank-neutral;
timing-dependent owner maps are not persisted or promised reproducible.
The metric covers instrumented gas operators, not particle work or all MPI
communication. This is functionality evidence, not a scaling result.

## Completed focused verification

Reused GNU Fortran 13.3.0 Debug/Ninja fixed tests-enabled and selected H2/O2
tests-disabled builds, `/usr/bin/gfortran`, Open MPI `/usr/bin/mpifort` and
`/usr/bin/mpiexec`. Compiler/MPI/backend configuration did not change.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_balance pelef_mpi_amr_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_balance_(2|4)$' --output-on-failure --output-log build/spray-amr-balance-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_mpi_amr_spray_restart.py --exe build/spray-amr-reactive-debug/pelef_mpi_amr_spray_3d --mpiexec /usr/bin/mpiexec --source /home/gtt13/projects/pelec --work /home/gtt13/projects/pelec/build/spray-amr-balance-checks/fixed --multiport --adaptive --rebalance
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_spray_3d_selected --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_mpi_amr_spray_restart.py --exe build/spray-amr-reactive-selected-debug/pelef_mpi_amr_spray_3d_selected --mpiexec /usr/bin/mpiexec --source /home/gtt13/projects/pelec --work /home/gtt13/projects/pelec/build/spray-amr-balance-checks/selected --selected --multiport --adaptive --rebalance
python3 -m py_compile tests/application/check_mpi_amr_spray_restart.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final deterministic balancing CTest **2/2 passed**, no skips: two ranks
**0.55 s**, four ranks **0.57 s**, **1.16 s** total. The fixture concentrates
all three patches/particles on rank zero, supplies known unequal positive
costs, verifies redistribution, exact gas/temperature values, particle
composition and finest-owner assignment. Receive capacity zero forces a late
particle migration failure and verifies owner/U/T/particle rollback. Zero
costs verify a successful no-op. The fixture also retains its remap prerequisites.

Both fixed and selected measured-time application checkers passed, including
adaptive two-port runs, two-rank save/four-rank continuation, gas/particle CSV
agreement at relative 2e-12 / absolute 1e-25, future injection events and the
existing context/truncation/overwrite rejection checks. The application checker
is also registered as `spray_mpi_amr_rebalance_application`; the completed run
here used the direct Python command above.

Initial fixed/selected builds failed because the extended namelist declaration
exceeded the Fortran line-length limit. The initial command batch continued:
CTest could not launch the unbuilt test, and stale application binaries rejected
the new namelist controls. Splitting the source line and requiring build success
before the retry resolved these failures; only the final runs above are passes.

Final provenance, Python syntax, project contract and whitespace checks passed.
Compiler warnings remain. No installed audit, timing reproducibility guarantee,
particle-cost model, independent energy/order oracle, scaling or physical
validation was performed. CLI fixtures disable gas chemistry/transport.

Logs: `build/spray-amr-balance-checks/` preserves failed build/test attempts and
final `build-retry.log`, `build-selected-retry.log`, `ctest.log`, `test-fixed.log`
and `test-selected.log`. Application artifacts are under its
`{fixed,selected}/restart_*` directories.
