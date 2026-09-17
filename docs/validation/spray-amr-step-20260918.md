# Periodic MPI AMR mixed-particle intervals — 2026-09-18

Implementation: `17fb7577a36dfe396c17caa15356f88aca8618ec`, based on
`7d5ecab`. `advance_mpi_amr_mixture_coupling` connects existing mixture transfer
and breakup kernels to owner-local composite sampling and gas feedback.
Conserved fields and temperature guesses share the sample exchange. A common
minimum substep includes the finest tree spacing as its particle CFL bound.
Particle or gas-feedback rejection halves the trial step, up to 30 attempts.
Feedback uses pre-step positions; trajectories wrap and survivors migrate at
the end of the interval. U/T, parcels, component inventories and optional ID
watermark publish together only after final migration succeeds.

The caller supplies common species/transport/liquid tables and options, globally
unique locally increasing particle IDs, and a static periodic tree for the
interval. ETAB/KH-RT representative breakup and KH children reuse existing
models; child transverse kinetic energy is debited from the gas. This library
does not yet orchestrate gas hydro/chemistry, regridding, injection or restart.

## Completed focused verification

GNU Fortran 13.3.0, Debug/Ninja, explicit `/usr/bin/gfortran`, Open MPI
`/usr/bin/mpifort` and `/usr/bin/mpiexec`. Reused the fixed tests-enabled and
selected tests-disabled build directories from the immediately preceding
[gas-feedback verification](spray-amr-mpi-gas-20260918.md); compiler, MPI and
backend configuration did not change. No installed prefix was used.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-mpi-gas-debug --target pelef_spray_tree_runtime --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-mpi-gas-debug --target test_spray_mpi_amr_step --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-mpi-gas-debug -R '^spray_mpi_amr_step_(2|4)$' --output-on-failure --output-log build/spray-amr-step-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-mpi-gas-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final **2/2 passed**, no skips: two ranks **0.75 s**, four ranks **0.68 s**,
**1.47 s** total. The three-level water/peroxide fixture starts all particles on
one rank, exercises coarse/fine crossings and periodic wrapping, verifies
finest-owner migration and active drag/heating/evaporation, and checks composite
gas plus liquid mass, momentum, energy and species conservation within 3e-12
of each inventory scale. A separate KH-child case checks ID advancement,
child survival and the same coupled conservation budget. Empty populations,
interval disagreement, late substep exhaustion and final receive-capacity
rejection cover collective behavior and whole-interval rollback.

The first run failed both tests because the test's 2 microsecond interval fit
one substep, so its exhaustion expectation was false. Changing the fixture to
20 microseconds forces multiple steps through the particle CFL bound, without
changing solver tolerances. The base checks then passed 2/2 in 1.37 s; adding
KH-child and empty-population coverage produced the final result above. Initial
Open MPI output also contained display-authorization messages; processes ran
and the reported test failures were the explicit exhaustion assertions.

The selected H2/O2 library compiled successfully; this is compile evidence,
not selected executable or installed-binary qualification. Provenance, Python
syntax, project contract and whitespace checks passed. Compiler warnings remain.
No Python source changed. No external solver parity, ETAB-specific new test,
allocation-failure injection, AMR convergence, scaling or physical spray
validation was performed. Physical particle boundaries, AMR wall films and EB
collisions remain open.

Logs: `build/spray-amr-step-checks/{build,build-test,build-selected,ctest-initial,ctest-base,ctest}.log`.
Executable: `build/spray-amr-mpi-gas-debug/tests/test_spray_mpi_amr_step`.
