# Coupled distributed EB remap — 2026-09-18

Implementation `357dd763778e03745fa53518b6dd3f6540ce08e1`, base `6aed2c7`.

`remap_mpi_eb_spray` wraps the existing conservative MPI EB gas remap and EB
particle migration in private state. Topology, owners, geometry, gas/temperature,
parcels and component masses publish only after both stages succeed. Gas uses
the existing shared-root/refinement-lattice and immutable-geometry contract.
The caller provides consistent gas geometry and original nodal values.

The particle stencil, migration and MPI sampling geometry checks now allow an
original nodal lattice that is a common positive integer refinement of the
current finest tree lattice. This retains the original collision/visibility
geometry after gas coarsening. Dimensions remain positive, individually bounded,
and at most one million cells in total. MPI sampling explicitly agrees nodal
array dimensions before broadcasting values, including ranks with no queries.
No geometry is reconstructed from aggregated cut-cell moments.

## Completed focused checks

GNU Fortran 13.3.0 `/usr/bin/gfortran`, Open MPI 4.1.6 `/usr/bin/mpifort` and
`/usr/bin/mpiexec`, Debug/Ninja. Existing fixed tests-enabled and selected H2/O2
tests-disabled configurations were reused unchanged.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_remap test_spray_mpi_eb_sample test_spray_mpi_eb_tree --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_(remap_[24]|sample_2|tree_2)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py tools/check_mechanism_provenance.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The initial focused CTest set passed **4/4**, total **2.01 s**: remap two/four
ranks 0.51/0.54 s, existing sampling 0.47 s and migration 0.45 s. After adding
an explicit mismatched-rank nodal-dimension rollback assertion, the final remap
rerun passed **2/2**, 0.51/0.52 s, total **1.07 s**. These are four distinct tests,
not six. Selected compilation, provenance, syntax, contract and whitespace
passed. No failures, skips or tolerance changes; existing compiler warnings
remain.

The fixture uses a 4 mm cube, 4³ root, an 8³ nodal plane y=1.25 mm, an initially
left refined child, and one mixed water/peroxide parcel on the fine owner.
Fine gas density/conserved state is multiplied by 1.1 so coarsening requires a
real weighted average. It moves refinement right, coarsens to root-only, then
refines left again, changing owners along the way. Composite gas inventories
remain within relative 3e-13 (floor 1e-30); particle identity, per-drop mass and
components are preserved. Each stage samples the gas using the original 8³
nodal field. On the root-only stage, one rank supplying a 4³ nodal field while
others supply 8³ rejects collectively and preserves prior samples.

A zero particle capacity forces rejection after private gas remap; owner map,
gas, temperature and particle/component state remain unchanged. This tests late
coupling failure, not only rejection before gas work.

Artifacts: `build/spray-mpi-eb-remap-checks/{build.log,build-final.log,build-selected.log,ctest.log,ctest-final.log}`.

This is a library remap connection. The EB spray application's regrid scheduling
and checkpoint-stage integration remain disabled pending connection. Selected
compilation is not selected runtime evidence. No convergence, performance,
external parity, full matrix or physical validation was run.
