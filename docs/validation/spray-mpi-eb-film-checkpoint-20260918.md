# MPI EB particle and film checkpoint — 2026-09-18

Implementation `37111e6`, base `cb56173`.

The checkpoint API accepts optional films, film component masses and capacity.
It retains V1 for particle-only calls and uses V2 when films are supplied.
The writer gathers private films through EB migration with all tree owners set
to rank zero, then publishes one immutable particle/film archive. The reader
stages all payloads, metadata and owner redistribution before committing either
population. The caller still binds gas, geometry, models and schedules through
its generation context.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, `/usr/bin/mpifort`, `/usr/bin/mpiexec`,
Debug/Ninja; reused fixed tests-enabled and selected H2O2 tests-disabled builds.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_checkpoint --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_particle_checkpoint$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tests/application/check_mpi_eb_particle_checkpoint.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **1/1 passed**, test time **3.39 s**, total **3.43 s**, no failures or skips.
The driver performs seven MPI launches: V1 write/read/truncation and V2
write/read/truncation/invalid-normal. Writers use two ranks and readers use four.
Each writer also checks immutable-path rejection; each reader first checks a
wrong generation context. Selected library compilation, provenance, Python
syntax, project contract and whitespace passed. Selected execution was not run.

The V2 fixture adds two mixed-component films, initially on the final writer
rank, on fine and exposed root portions of phi=y-0.00125. Restore places them
on their distinct current owners. Position, normal, mass, temperature, diameter,
height, material and component masses are checked exactly. Particle lineage,
deformation and KHRT state plus time, step, event, shared ID watermark and
injection inventory retain the existing exact checks.

For the late-failure case the driver reverses one stored film normal while
retaining a unit-length normal and otherwise valid payload. Particle reading
and redistribution proceed privately, but EB film migration rejects the surface
orientation. Empty caller particle/film populations and metadata sentinels
remain unchanged. Truncated files and wrong generation contexts also roll back.

Artifacts: `build/spray-mpi-eb-film-checkpoint-checks/{build.log,ctest.log,build-selected.log,provenance.log,contract.log}`;
per-launch logs and archives are below
`build/spray-amr-reactive-debug/tests/spray_mpi_eb_particle_checkpoint/`.

The capacity also bounds the root-gathered film archive. Films do not acquire
stable IDs or globally canonical ordering. Gas state remains a separate payload.
MPI application film input/output and generation-manifest integration remain
open. No broad matrix, installed audit, external parity, scaling or physical
validation was run.
