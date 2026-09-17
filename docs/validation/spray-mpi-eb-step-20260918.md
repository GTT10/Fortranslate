# Distributed EB mixed-particle interval — 2026-09-18

Implementation: `9e99cb9e3780ad7a3e470dae607a7bae0c21e001`, based on `5456221`.

Paired optional geometry/nodal fields select EB behavior in the existing MPI
AMR particle interval. Geometry mode agrees collectively, then EB sampling,
multicomponent source evolution, representative-drop breakup, source-first
specular drift and EB gas deposition operate on private state. Accepted local
wall exchanges accumulate and are globally summed once. Final EB owner migration
precedes publication of gas, temperature, parcels, component masses and ledger.
Failed trials halve the common substep; failed intervals preserve caller state
and return zero accepted substeps. Regular AMR remains the default path.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja; reused unchanged fixed
tests-enabled and selected H2/O2 tests-disabled configurations.
Compiler `/usr/bin/gfortran`, MPI wrapper `/usr/bin/mpifort`, launcher
`/usr/bin/mpiexec`.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_step test_spray_mpi_amr_step --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_mpi_eb_step_[24]|spray_mpi_amr_step_2)$' --output-on-failure --output-log build/spray-mpi-eb-step-checks/ctest-final.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tools/check_project_contract.py tools/check_mechanism_provenance.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **3/3 passed**, total **1.83 s**: EB two ranks **0.53 s**, EB four
ranks **0.55 s**, regular two-rank interval regression **0.72 s**. Selected
compilation, provenance, syntax, contract and whitespace checks passed. No
failures or skips; existing compiler warnings remain. An initial successful EB
fixture placed the particle on rank zero; the final fixture deliberately starts
it on the fine owner and requires migration to the coarse owner. These are
revisions of the same tests, not additional distinct cases.

The final fixture has a 4 mm cube, 1 mm root cells, one ratio-two child and EB
plane y=1.25 mm. Gas is N2 at 1000 K with nonzero momentum. One 350 K mixed
water/peroxide parcel starts 1 nm from the wall near the coarse/fine interface,
with velocity (20,-20,5) m/s, mass 1e-12 kg and multiplicity 100. A 20 microsecond
interval requires multiple source substeps. Assertions require heating,
evaporation, drag, reflection, positive wall impulse, crossing to an exposed
coarse cell and final migration from rank N-1 to rank zero. All composite
gas-plus-liquid components close against the wall ledger at relative 3e-12 with
1e-30 scale floor. A one-substep limit rejects after partial private progress and
preserves incoming gas, temperature, particle population/position, components
and ledger. Breakup is disabled in this new fixture.

Artifacts: `build/spray-mpi-eb-step-checks/` contains `build-final.log`,
`build-selected.log`, `ctest-final.log` and the initial fixture logs.

EB mode rejects explicit exterior particle face settings, wall-film arguments
and KH children. It inherits root-contained drift and first-order source/drift
splitting. The gas receives sources but its hydro/chemistry/transport operators
are not advanced by this API. Distributed full gas/spray composition, executable
injection/restart and physical wall films remain open. No convergence,
production scaling, external parity or physical validation was performed.
