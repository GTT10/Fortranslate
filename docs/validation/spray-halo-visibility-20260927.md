# EB visibility in particle halo sampling — 2026-09-27

Implementation `1053688d`, based on `befa7952`, local `main`, version 0.246.0.
PeleC reference `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance preflight SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Change

The existing composite EB stencil's periodic segment-splitting visibility
routine is exposed as `eb_particle_centroid_blocked` and reused by particle
halo interpolation. Optional paired `nodal_phi` and `centroids` arguments
require a fluid-fraction halo. Centroids contain absolute physical coordinates
on the halo layout. Positive-weight fluid donors blocked by the nodal wall
are removed, followed by weight normalization. Solid queries and zero remaining
support reject collectively without changing samples.

The sampler agrees on option presence and nodal dimensions, validates lattice
compatibility/finiteness, compares the replicated nodal field, and checks
periodic endpoint equality. The stage-time gas wrapper forwards the pair and
requires explicit target-volume fractions in this mode. Existing gas halo
validation retains exact interior metric agreement. Centroid shapes and
requested values are checked; their provenance from the corresponding extended
geometry remains a caller contract. Omitting the pair retains fraction-only
sampling. Source transfer and recursive particle application stepping remain
unconnected.

## Completed verification

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6 through `/usr/bin/mpifort` and `/usr/bin/mpiexec`.

Six CTests passed, 3.23 s final run:
`spray_mpi_amr_sampling_2`, `_4`, `spray_mpi_eb_sample_2`, `_4`,
`spray_mpi_eb_periodic_sample_2`, `_4`.

A manufactured nodal slab separates two otherwise enabled interpolation donors
with values 10 and 100. Sampling from the left yields only 10 (second component
-20), within 1e-13. Moving a query into the solid or changing nodal geometry on
one rank rejects and retains the prior samples. This is an isolated visibility
fixture; its unit donor fractions do not represent a physical slab geometry.
The temporal EB fixture now supplies explicit planar target fractions and
centroids and enables visibility for solver-generated halos at parent fractions
0, 0.25 and 1. Both ordinary and periodic configurations pass. Existing capped,
composite and halo sampling checks also execute.

Initial compilation rejected an overlong declaration; it was wrapped. The
first temporal fixture recomputed interior planar volume fractions and was
rejected by the gas halo's exact metric identity check. Interior fractions now
come directly from the existing geometry; tolerances were not loosened. Initial
failure logs are retained. Existing compiler warnings and X authorization
messages remain. Provenance, final target builds, Python syntax, project contract
and whitespace checks passed.

No full matrix, selected-mechanism build, installed audit, recursive application
run, arbitrary-wall physical validation or upstream trajectory parity was run.
Periodic wall-crossing occlusion itself was not independently isolated in this
increment; the shared periodic intersection path has existing composite tests.

## Reproduction

Artifacts: `build/spray-halo-visibility-checks/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_sampling test_spray_mpi_eb_sample --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_(amr_sampling|eb_(periodic_)?sample)_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
