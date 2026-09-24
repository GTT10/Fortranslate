# MPI LES donor-plane sampling — 2026-09-25

Implementation `cb8a898`, based on `e247e3f`.
`sample_mpi_les_x_planes_3d` fetches supplied one-based global x coordinates
from the corresponding owners of an uneven x-slab decomposition. Each rank
may request 0..256 planes, including duplicates. Integer requests are routed
with Alltoallv; packed field planes return in a second Alltoallv and are restored
to the caller's request order. No global field array is assembled. The caller
owns physical-coordinate mapping and wall transforms; this API transfers raw
in-domain donors for the pending regular MPI Dynamic boundary integration.

Grid/component configuration, local shapes, finite fields, coordinate bounds,
MPI count products and allocations are checked collectively. The output is
allocatable and replaced only after exchange completes successfully. Zero
requests are supported on a rank that still owns a nonempty slab. Empty spatial
slabs are not supported. Repeated requests are preserved rather than deduplicated.

## Completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpiexec`), Ninja Debug. Fixed unit
and selected H2O2 spray runtime built. Selected execution was not run.
CTest `spray_les_mpi_filter_2` and `_4`: **2/2 passed**, 0.47 and 0.49 seconds,
total 1.00 second. No failures, skips or tolerance changes.

A 5x3x2 field is split into uneven slabs. Twenty-six periodic/reflected donor
coordinates span multiple slab widths and match serial planes exactly.
Rank-dependent request counts include zero and repeated last-owner coordinates.
A coordinate invalid on one rank rejects collectively and preserves previous
plane values. Existing box/Gaussian/identity filter parity and rejection checks
also pass. These checks verify transport of already-mapped coordinates, not
physical reflection of momentum or Dynamic boundary flux integration.
Mechanism provenance, Python syntax of the contract checker, project contract
and whitespace passed. This is MPI communication evidence, not a spray run,
physical validation or external solver parity.

## Commands and artifacts

Logs: `build/mpi-les-plane-checks/`. Launcher overrides are host-local.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_mpi_filter --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_mpi_filter_(2|4)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
