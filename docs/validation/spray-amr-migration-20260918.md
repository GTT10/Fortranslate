# MPI spray migration on AMR trees — 2026-09-18

Implementation: `07d43b4876967af131c071041dca8fbd68300724`, based on
`bd13c51`. `migrate_amr_tree_spray_parcels` validates common topology, owner map,
physical root bounds/spacing, periodic flags and component mode. It partitions
private parcel/component arrays using the finest-cell owner, exchanges the
rank bins and verifies received ownership before publishing both arrays.

The Cartesian migrator and tree migrator now share `exchange_spray_parcel_bins`:
count exchange, bounded payload allocation, integer/real Alltoallv, ID sorting,
received-state/component validation and collective commit. Parcel IDs, parent
IDs and material indices retain integer transport; all fourteen stored real
values and optional component columns retain their prior transport layout.
The existing Cartesian gather option continues through the same helper.

Local capacity can differ by rank, from zero to one million parcels. Component
presence/count must agree. Metadata, capacity and received validation failures
preserve original application arrays on every rank. IDs must be globally unique
by caller contract; the receive check rejects duplicates on one destination,
not duplicates deliberately assigned to different destinations.

This is migration of an accepted particle population. It does not yet connect
tree gas sampling/deposition, coupled particle stepping, checkpoint integration
or EB impacts. Periodic lookup does not change stored positions.

## Completed focused verification

Two clean GNU Fortran 13.3.0 Debug/Ninja configurations with `/usr/bin/gfortran`,
Open MPI `/usr/bin/mpifort` and `/usr/bin/mpiexec`: tests-enabled fixed runtime
and tests-disabled selected H2/O2 library build. No installed prefix is used.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/spray-amr-migration-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_SPRAY=ON -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=ON
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-migration-debug --target test_spray_mpi_amr_tree test_spray_mpi --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-migration-debug -R '^(spray_mpi_amr_tree_(2|4)|spray_mpi_migration_2)$' --output-on-failure --output-log build/spray-amr-migration-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/spray-amr-migration-selected-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_SPRAY=ON -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=OFF -DPELEF_MECHANISM_BUNDLE=mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-migration-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**3/3 passed**, no failures/skips, **9.64 s** total:

- AMR tree migration, two ranks: **0.46 s**.
- Existing Cartesian MPI migration/coupling regression, two ranks: **8.65 s**.
- AMR tree migration, four ranks: **0.50 s**.

The new three-level fixtures verify exact complete parcel/component payloads,
receive ID ordering, empty destinations/populations, periodic finest ownership,
particle crossings, changed gas ownership, coarse-only regrid migration and
the payload without component columns. Capacity overflow, inconsistent owner
maps/root bounds and duplicate received IDs reject with rollback. The existing
Cartesian regression also exercises its gather/checkpoint and coupled consumers
after exchange extraction.

The selected library compiled successfully; its Ninja dependency edge contains
selected libraries without `pelef_core`. This is compile/dependency evidence,
not selected executable or installed-binary evidence. Provenance, Python syntax,
project contract and whitespace checks passed. No failed build/test attempt;
compiler warnings remain. No Python source changed. No memory-failure injection,
MPI fault recovery, scaling, coupled AMR gas/particle or physical validation is
claimed.

Logs: `build/spray-amr-migration-checks/{configure,build,ctest,configure-selected,build-selected}.log`.
Test executables: `build/spray-amr-migration-debug/tests/{test_spray_mpi_amr_tree,test_spray_mpi}`.
