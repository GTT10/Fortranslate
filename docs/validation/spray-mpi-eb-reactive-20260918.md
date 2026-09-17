# Distributed EB gas/spray composition — 2026-09-18

Implementation `4a25e3845fa8ef1b1e07bbbe3bceaab8873a1f4f`, base `1c7999d`.

The new wrapper composes particle half-steps around the existing MPI EB full
gas operator. Owner-local gas/temperature, particles and component masses are
private until both halves and gas evolution succeed. Exterior, wall, chemistry
and limiter diagnostics publish at the same point. Particle wall impulses are
added to the gas wall ledger. Failure returns zero accepted particle steps.
Selected runtime explicitly links the selected MPI EB tree runtime.

## Completed focused checks

GNU Fortran 13.3.0 (`/usr/bin/gfortran`), Open MPI 4.1.6
(`/usr/bin/mpifort`, `/usr/bin/mpiexec`), Debug/Ninja. Existing fixed tests-enabled
and selected H2/O2 tests-disabled build configurations were reused unchanged.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_reactive --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_reactive_(2|4)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py tools/check_mechanism_provenance.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **2/2 passed**, two ranks **1.01 s**, four ranks **1.02 s**, total
**2.06 s**. Selected compilation, provenance, syntax, contract and whitespace
checks passed. No failures or skips; existing compiler warnings remain.
An initial successful run preceded explicit initialization of private diagnostics
and stronger parcel rollback assertions; it is not counted as extra coverage.
Logs: `build/spray-mpi-eb-reactive-checks/{build-final.log,ctest-final.log,build-selected-final.log}`.

The fixture is a 4 mm cube, ratio-two left child, plane y=1.25 mm, N2 at 1000 K,
and one water/peroxide mixed parcel on the fine owner. A 10 ns step enables
implicit chemistry and all four molecular transport processes. It requires
nonzero chemical species exchange, particle reflection, heating, evaporation
and drag. Composite gas-plus-particle inventories close against exterior,
wall and chemical ledgers within relative 3e-12 with a 1e-30 scale floor.
An invalid gas solver rejects after the first particle half-step and preserves
owner-local gas/temperature, particle population, mass, position, velocity,
temperature, component masses and all diagnostics exactly. Empty intermediate
ranks participate in the four-rank case.

This focused test uses common gas steps, no child creation and no breakup.
Optional gas subcycling is forwarded but not independently exercised here.
Selected compilation is not a selected runtime test. This wrapper inherits
root-contained specular drift, first-order source/drift splitting and the EB
composite stencil approximation. Distributed executable injection/restart,
liquid wall films, convergence, scaling, external parity and physical validation
remain open. It does not establish complete PeleC parity.
