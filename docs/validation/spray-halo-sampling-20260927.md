# Particle sampling from gas halos — 2026-09-27

Implementation `ecf6ce43`, based on `21ddab24`, local `main`, version 0.246.0.
PeleC reference `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance preflight SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Interface and boundary

`sample_mpi_amr_particle_halo` interpolates caller-filled full gas halos using
assigned target-node indices. Query nodes must be locally owned. The routine
agrees on topology, ownership, component count, width, periodic axes and
physical bounds before sampling. It supports native spatial lower bounds and
packed storage of equal extents. Trilinear weights use target-resolution cell
centers and the periodic image nearest the patch center. Optional matching
fluid-fraction halos multiply and renormalize weights. Zero support, missing
halo extent, invalid requested values or configuration disagreement preserves
all sample outputs collectively.

The existing `fill_mpi_eb_tree_halo_3d` produces this layout and accepts ancestor
endpoints/fractions. Its volume-fraction output requests complete corner donors;
wider EB halos require explicit target metrics. The new sampler consumes such
arrays but does not itself call halo filling, interpolate endpoints in time or
check EB centroid visibility. These connections and visibility protection are
still required before arbitrary cut-cell particle coupling. Recursive spray
source transfer and application stepping remain open.

## Completed verification

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6 through `/usr/bin/mpifort` and `/usr/bin/mpiexec`.

`spray_mpi_amr_sampling_2` and `_4` passed (2/2, 1.05 s total).
The new fixture constructs two-layer halos with an affine two-component field
`[f,-2*f]`, `f=2*x+3*y-z+7*0.25`. Queries just beyond a fine patch and across
a periodic boundary match analytic values within 1e-13. A uniform 0.5 fluid
fraction leaves the normalized result unchanged. Native and packed halo lower
bounds give exact matching outputs. Zero fluid support, rank-disagreeing halo
width and a query outside the available halo reject and retain prior samples.
Existing composite/capped sampling assertions also execute.

The time-like term is manufactured input; this is not a runtime check of
ancestor time interpolation. Actual solver-generated halo consumption,
nonuniform cut-cell fractions, EB line-of-sight and temporal convergence were
not independently tested in this increment.

Provenance, target build, Python syntax, project contract and whitespace checks
passed. No build/test failures occurred. Existing compiler warnings remain.
No full matrix, selected-mechanism build, installed audit, recursive application
run or physical validation was performed.

## Reproduction

Artifacts: `build/spray-halo-sampling-checks/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_sampling --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_sampling_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
