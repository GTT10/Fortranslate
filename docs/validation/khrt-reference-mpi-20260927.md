# Variable-count KHRT release APIs — 2026-09-27

Implementation `9a0130e5`, based on `f5ad03bf`, local main, version 0.246.0.

## Change

`advance_mixture_khrt_shedding` and `advance_mpi_khrt_shedding` accept paired
optional arguments `reference_factor` and allocatable `gas_momentum(3,np)`.
When present, each release uses the reference child-count/azimuth constructor
from [the preceding kernel increment](khrt-reference-children-20260927.md).
The default four-child path retains its existing API and behavior. A missing
partner argument is rejected so callers cannot silently discard momentum.

Counts are calculated before allocation and ID reservation. Capacity rejection
preserves parents, component masses, ID watermark and both feedback outputs.
MPI agrees on mode and exponent, including idle ranks. The existing parent-ID
ordering handles variable group sizes and assigns rank-independent child IDs.
Momentum and energy increments are returned at original parent sites; callers
must deposit them in the same gas/particle transaction.

`khrt_reference_uniforms` generates stateless xorshift64 samples keyed by parent
ID and the pre-release global ID watermark. Those identifiers already persist
through checkpoint/restart. Bit shifts/XOR avoid signed arithmetic overflow;
53-bit fractions lie in [0,1). A rejected candidate advances no random state.
This is a PeleF random sequence, not reproduction of AMReX's random stream.

## Independently completed checks

Existing fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`,
GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpifort`, `/usr/bin/mpiexec`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi test_spray_core --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_(mpi_migration_[24]|khrt_reference_children)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final **3/3 passed**, no failures/skips: MPI 2 ranks 10.73 s, MPI 4 ranks
11.38 s, constructor unit 0.00 s; total 22.17 s. Build, provenance, Python
syntax, contract and whitespace checks passed. Existing compiler warnings
remain. Artifacts: `build/khrt-reference-mpi-checks/`, including preserved
`LastTest.log`. Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

The added MPI fixture distributes parents in reverse ID order with an idle
rank. Exponent 0.8 produces 27 children for one parent on 2 ranks, and 105
children across three parents on 4 ranks. Every child ID, velocity, mass,
multiplicity and component mass agrees exactly with a serial execution of the
same complete parent set. Per-parent momentum and energy feedback also agree
exactly. The existing 2e-12 inventory criterion passes after adding gas
compensation. Rank disagreement and capacity failure preserve caller outputs;
successful generation after those rejections still matches the serial result.
Existing migration/coupling checks in the same executables also pass.

## Remaining work

Gas deposition, coupled stepping, application controls and restart identity
still need to select and forward the reference policy. Current applications
continue to use four symmetric children. This check does not exercise a
changed-rank checkpoint, selected mechanisms, statistical random-distribution
quality, full upstream solver parity or physical validation.
