# MPI particle migration — 2026-09-16

Tested parent: `5ca6320` plus this change. Source/test hashes, configure/build/test
logs and resolved libraries are in `build/spray-mpi-migration-checks/`.

The new mpi_f08 API partitions outgoing particles, exchanges counts and separate
integer/real payloads, sorts received IDs, validates candidates and commits collectively.
Empty ranks and optional component columns are supported. Geometry/component modes
must agree across ranks; local input/capacity/received-record failures reject without
changing caller arrays. Process or MPI transport failure recovery is not promised.
IDs must be globally unique on entry; same-destination duplicate IDs reject.
Distributed gas feedback, film migration and application restart remain open.

## Focused verification

GNU Fortran 13.3.0 Debug with runtime bounds/FPE checks, Open MPI 4.1.6:

```sh
cmake -S . -B build/spray-mpi-migration-gnu-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec.openmpi -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_SPRAY=ON -DBUILD_TESTING=ON -DMPI_mpi_LIBRARY=/usr/lib/x86_64-linux-gnu/libmpi.so -DMPI_hwloc_LIBRARY=/usr/lib/x86_64-linux-gnu/libhwloc.so
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 8
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Two MPI tests (2 and 4 ranks) passed in 0.81 s total, no test failures/skips. Checks
cover gather then redistribution through empty ranks, exact ID/KH/component/velocity
payloads, received ID order, collective receive-capacity rollback and one-rank invalid
component rejection preserving other ranks.

The first clean configuration selected an incompatible Intel MPI wrapper and failed
the project's mpi_f08 preflight. A second clean GNU/Open MPI configuration initially
found Intel libmpi through the environment; explicit MPI library paths removed that
mix before tests. Resolved MPI libraries are Open MPI .so.40. Runtime hwloc remains
Intel's compatible libhwloc.so.15 via the host environment, recorded in the library
log. Provenance, contract and whitespace checks passed; no Python source changed.
No full matrix, distributed spray application or physical/scaling qualification ran.
