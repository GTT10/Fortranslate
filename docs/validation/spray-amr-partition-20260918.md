# Spray ownership on regular 3D AMR trees — 2026-09-18

Implementation: `0aa057383e46cf566ef3bbed86dc6cdd95efaab8`, based on
`08fbe44`. The new `spray_amr_tree_partition_mod` supplies batch finest-cell
lookup and stable rank-grouped parcel/component ordering. It uses existing
validated tree integer origins and level domains. Half-open cells and patches
give unambiguous peer/coarse-fine face ownership. Periodic axes wrap for lookup
without changing stored positions; nonperiodic exterior points fail.

The partitioner takes the current node owner map, so the same API can reassign
after coarsening, refinement or ownership changes. IDs must be strictly increasing
on input; sorting into rank bins preserves input order within each bin. Full
derived-type assignment preserves lineage, material, deformation and KH/RT
history. Optional component columns move with their parcels. Counts, offsets,
locations, parcels and components publish only after validation succeeds.

The separate opt-in `pelef_spray_tree_runtime` and
`pelef_selected_spray_tree_runtime` targets link their corresponding spray/tree
libraries. Existing spray executables do not acquire this dependency. No new
MPI migration, sampling, deposition, coupled AMR stepping or EB collision path
is claimed. Coordinate/spacing and integer-domain bounds are documented in
`docs/spray.md`.

## Completed focused verification

GNU Fortran 13.3.0, Debug/Ninja, `/usr/bin/gfortran`. Two clean configurations:
tests-enabled fixed runtime and tests-disabled selected H2/O2 library build.
No installed prefix is involved.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/spray-amr-partition-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DPELEF_ENABLE_SPRAY=ON -DPELEF_ENABLE_TESTS=ON
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-partition-debug --target test_spray_amr_tree_partition --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-partition-debug -R '^spray_amr_tree_partition$' --output-on-failure --output-log build/spray-amr-partition-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/spray-amr-partition-selected-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DPELEF_ENABLE_SPRAY=ON -DPELEF_ENABLE_TESTS=OFF -DPELEF_MECHANISM_BUNDLE=mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-partition-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 passed**, no failures/skips, reported test time **0.00 s**, total **0.03 s**.
The three-level fixture checks known finest nodes/local indices at nested and
peer boundaries, positive/negative periodic images, nonperiodic upper/exterior
rejection, stable rank bins with empty owners, exact preservation of every
stored parcel field and component association, invalid inventory/owner and NaN
rollback, duplicate-ID rejection, coarse-only reassignment and empty populations.

The selected library compiled successfully. Its Ninja dependency edge contains
only selected runtime libraries, without `pelef_core`; this is compile/dependency
evidence, not a selected executable or installed-binary audit. Provenance,
Python syntax, project contract and whitespace checks passed. No build/test
attempt failed; compiler warnings remain. No Python source changed. This is an
ownership/partition unit check, not coupled AMR particle conservation, MPI
migration, scaling or physical validation.

Logs: `build/spray-amr-partition-checks/{configure,build,ctest,configure-selected,build-selected}.log`.
Executable: `build/spray-amr-partition-debug/tests/test_spray_amr_tree_partition`.
