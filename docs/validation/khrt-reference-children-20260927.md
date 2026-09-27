# Reference KHRT child construction — 2026-09-27

Implementation `59f051b3e76d3318bb64962e4ee3cd3c7dc691c5`, based on `dd36ad97`,
local main, version 0.246.0.

## Reference and implementation

The local `build/breakup-reference/SpraySB.cpp` matches the SHA256 recorded in
`references/spray_khrt.json`:
`e540abf90a3e5acd9b91d6f7994c9373ad505bb47f512fa8c5b7c3e536f32382`.
It belongs to PelePhysics `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`.
Lines 133–136 select `N_d=max(1,int(N/N**factor))` and then `N_s=N/N_d`.
Lines 139–150 construct the transverse basis; lines 175–183 use a uniform
azimuth, add transverse velocity and drift from the saved release location.

`khrt_reference_child_count` implements the count formula with capacity and
integer-conversion checks. `construct_reference_khrt_children` accepts one
caller-supplied uniform sample in [0,1) per child. It implements the 3D basis,
velocity, equal multiplicities and initial KH count, retaining temperature,
material and lineage. It returns children at the release position; the caller
must apply the existing drift exactly once. The exponent must be in [0,1].

Random directions generally have nonzero net transverse momentum. To preserve
PeleF's conservation contract, the constructor additionally returns opposite
gas momentum and total-energy changes. These compensation terms are a PeleF
policy, not an upstream source-parity claim. Outputs commit together only after
validation and allocation succeed. No random generator or global ID counter is
modified internally, allowing caller-managed restart and MPI transactions.

## Independently completed checks

Existing Ninja Debug fixed-mechanism build `build/spray-amr-reactive-debug`,
GNU Fortran 13.3.0. The unit executable uses bounds/runtime checks and traps
invalid, zero-divide and overflow exceptions. No MPI execution was needed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_core --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_(khrt|khrt_reference_children)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final **2/2 passed**, no failures/skips, total 0.06 s. Provenance, build, Python
syntax, project contract and whitespace checks passed. Artifacts are in
`build/khrt-reference-children-checks/`, including the preserved `LastTest.log`.
The mechanism source SHA256 remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

The new focused check covers exponents 0, 0.5 and 1; truncation, subunit count,
capacity/integer overflow rejection; three nonsymmetric azimuth samples;
parallel-axis basis fallback; mass/momentum/energy conservation with gas
compensation; lineage and initial weight; invalid samples and ID-overflow
rollback. The existing KHRT test also passes. The first conservation assertion
failed after subtracting large unchanged thermal energies. The final test
checks unchanged material/temperature separately and directly sums the kinetic
energy change, retaining the original 1e-17 residual threshold. Initial failure
logs are retained. Existing compiler warnings remain.

## Remaining integration

Existing application paths still use four symmetric children. Connecting the
new constructor requires dynamic child counts in serial/MPI allocation,
rank-neutral random samples, vector gas feedback, configuration and restart
identity. Those are required follow-on work; this record does not claim an
application change, stochastic distribution qualification, full reference
solver parity or physical validation.
