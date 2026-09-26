# Conservative coarse routing of particle sources — 2026-09-27

Implementation `55a7de3b`, based on `41130eaf`, local `main`, version 0.246.0.
PeleC reference `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance preflight SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Transfer

Optional `route_to_coarse` in `sum_mpi_amr_particle_source_halo` defaults false.
True routes in-domain contributions without a same-level receiver to the
finest covering ancestor. Integer level coordinates choose the destination.
Source densities are scaled by source/destination geometric cell-volume ratio.
Paired optional `eb_geometry` and `halo_volume_fraction` supply destination and
source fluid fractions, giving a fluid-volume ratio instead. Replicated EB
hierarchy validation and owner-local fraction shape/range checks precede routing.
Fractions must come from the source halo's target-resolution geometry; their
provenance is a caller contract.

Routed nonzero contributions with zero source/destination fluid volume, overflow
or zero scaling reject without publication. State and remainder remain atomic.
Nonperiodic domain-exterior entries stay in the remainder. Default false retains
unmatched coarse/fine entries. No parent restriction, gas EOS update, population
de-duplication or temporal scheduling is performed. This is conservative source
inventory transfer, not a claim of PeleC coarse/fine scheduling parity. Recursive
particle/film clocks and source synchronization still need driver integration.

## Completed verification

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6 through `/usr/bin/mpifort` and `/usr/bin/mpiexec`.

`spray_mpi_amr_deposition_2` and `_4` passed (2/2, 1.06 s total).
The added two-level fixture has a partial fine patch on a different rank from
the root. A fine-halo source density 12*[1,-2] lies beyond that patch. Default
routing leaves it in the remainder; geometric coarse routing gives 1.5*[1,-2]
with volume ratio 1/8. With source fluid fraction 0.5 and destination fraction
0.75, coarse routing gives [1,-2] within 1e-13. Together with a same-level source
2*[1,-2] in a full fine cell, the routed fluid-volume inventory is [1,-2]
within 1e-13. The routed remainder clears while a physical-exterior entry stays.
A nonzero contribution assigned zero source fluid volume rejects and retains
both outputs. Existing same-level periodic/peer routing, capacity rollback,
composite and capped deposition tests also execute.

Provenance, target build, Python syntax, project contract and whitespace checks
passed. No build/test failures occurred. Existing compiler warnings remain.
No full matrix, selected-mechanism build, installed audit, recursive application
run, temporal-order check or physical validation was performed.

## Reproduction

Artifacts: `build/spray-coarse-source-checks/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_deposition --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_deposition_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
