# MPI EB particle domain-wall reflection — 2026-09-26

Implementation 1b8c539, parent 75d0f96, local main, version 0.246.0.
The MPI EB mixed-particle interval now accepts explicit slip/no-slip outer
faces and forwards them to the existing level-set drift. Both use specular
particle reflection. Outer-wall momentum change is returned through the
existing exterior ledger. The reacting wrapper now constructs particle face
labels for wall configurations as well as periodic ones; gas inflow still
maps to particle outflow. Application input remains guarded until coupled
runtime checks are complete. Films and KH children with explicit EB particle
faces remain rejected by the existing interval validation.

## Completed check

GNU Fortran 13.3.0, Open MPI 4.1.6, Python 3.12.3, Ninja Debug. Existing fixed
tests-enabled build/spray-amr-reactive-debug, no backend changes.
Four-rank wall particle test passed **1/1**, 0.55 s (CTest total 0.60 s), no
failures or skips. Two gas patches leave two ranks empty. A parcel near x-high
reflects back into the coarse patch while also reflecting off the embedded
wall, heating and evaporating. Exterior mass/species/energy change is exactly
zero; x impulse is negative. Gas/particle change equals embedded and exterior
ledgers within the existing 3e-12 relative inventory bound. Existing forced
substep exhaustion verifies whole-interval rollback before the successful run.
The two-rank test entry was added but not run. Reacting wrapper forwarding was
compiled, not independently runtime-qualified in this increment.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_step --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_wall_step_4$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, project contract and whitespace passed. Compiler warnings
remain. No tolerance or golden-data changes.
Artifacts: build/mpi-eb-wall-particle-checks/build.log and LastTest.log.
Mechanism SHA: 0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4.
Pinned PeleC: bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8.
Pinned PelePhysics: 6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b.

## Limits

Coupled reacting runtime, wall application input/restart, films and KH children
remain separate work. No full matrix, installed audit, external solver parity
or physical spray-wall validation ran. Complete PeleC parity remains unproven.
