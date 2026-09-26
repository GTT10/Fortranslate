# Same-level particle source halo summation — 2026-09-27

Implementation `34cb41dd`, based on `5f050e2d`, local `main`, version 0.246.0.
PeleC reference `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance preflight SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Transfer behavior

`sum_mpi_amr_particle_source_halo` routes nonzero valid/halo contributions to
matching valid cells on a requested level. Integer coordinates handle patch
boundaries and periodic wrapping. MPI exchanges target-cell records and source
vectors; additions use increasing source-rank order with overflow rejection.
Other levels are unchanged and no parent restriction runs.

Source and remainder outputs publish together. Routed entries become zero in
the remainder; unmatched coarse/fine or physical-domain exterior entries stay
in their original halo slots. Input halos remain read-only. Caller source
values must be additive contributions with consistent normalization/EB volume
metrics for every image of a physical cell. Synchronized duplicate ghost data
must not be supplied as independent contributions.

Pinned PelePhysics `SprayParticles.H::transferSource` calls `SumBoundary` in
the PeleC path before adding temporary sources to active sources. This change
supplies same-level additive exchange. It does not implement unmatched-source
coarse/fine reflux, wall transformations, volume scaling, virtual/ghost
population accounting or recursive particle stepping. Callers must handle the
remainder and consume each halo once before claiming completed transfer.

## Completed checks

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6 through `/usr/bin/mpifort` and `/usr/bin/mpiexec`.

`spray_mpi_amr_deposition_2` and `_4` passed (2/2, 1.04 s total).
A two-peer fine level assigns left and right patches to different ranks.
Five input contributions use vector multiples of [1,-2]. Interior, opposite
peer halos and the periodic domain seam produce exact destination values:
left cell one receives 5*[1,-2], left cell four 3*[1,-2], and right cell one
2*[1,-2]. A nonperiodic exterior contribution remains as 5*[1,-2]. Thus source
plus remainder exactly retains all 15*[1,-2] of input. The root's sentinel
field remains unchanged. Zero receive capacity retains source and remainder;
rank-disagreeing widths reject. Existing composite/capped deposition tests run.

The coarse/fine-gap remainder case was not separately exercised; the tested
unmatched path uses a physical-domain exterior. No EB metric scaling is tested
because this routine performs none. Provenance, target build, Python syntax,
project contract and whitespace checks passed. No build/test failures occurred.
Existing compiler warnings remain. No full matrix, selected-mechanism build,
installed audit, recursive application run or physical validation was performed.

## Reproduction

Artifacts: `build/spray-source-halo-checks/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_deposition --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_deposition_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
