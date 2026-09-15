# Normal/lognormal diameter laws — 2026-09-16

## Source and implementation

Tested parent: `7b226d487a814108e5f39e72c6affd911045efbe` plus this commit's
source/tests. Exact tested file hashes and logs are retained in
`build/normal-size-checks/`.

The pinned PelePhysics distribution source and license are recorded in
`references/spray_size_distribution.json`. Normal and LogNormal now connect
to the Fortran injector, including single, shared multi-port and independent
port schedules. Arithmetic mean/std conversion follows the pinned source.
A small-ratio expansion avoids cancellation when deriving lognormal width.

The new inputs are `diameter_distribution` and `diameter_std`. Existing bounds
truncate the CDF; deterministic midpoint quantiles are inverted by bisection.
Positive-tail normal intervals use survival probabilities. Numerically
unresolved intervals are rejected. The existing mass/number weighting selects
normalization, with prescribed dose retained. Global standard deviation applies
to each port's own mean diameter. Distribution/std enter restart identity;
default Weibull/monodisperse contexts remain unchanged.

## Focused checks

GNU Fortran 13.3.0, existing Debug spray build, MPI disabled:

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d test_spray_core --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(size_distribution|application_cli)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Two CTest entries passed in 31.93 s, including 25 CLI tests; no failures or
skips. Build, provenance, Python syntax and contract checks passed. Compiler
warnings include real comparisons and existing bounds-analysis diagnostics.

For each new law, the CLI compares all eight batch diameters with independent
Python `statistics.NormalDist.inv_cdf` values (relative tolerance 2e-13),
checks equal number weights and conservation, then compares byte-identical
final gas/parcel CSVs after restart. Single, shared two-port and independent
schedule paths are covered. Changed std on restart, unknown distribution,
negative std and conflicting/unused parameters are rejected. Existing
Weibull mass/number tests also pass.

## Limits

Deterministic bounded quadrature and exact-dose normalization differ from
upstream random sampling. Truncation changes the actual population moments.
The source laws are implemented; this does not establish stochastic trajectory
parity or physical spray validation. Chi-squared sizes, physical child
formation/dispersion, mixture droplets and particle MPI/AMR/EB remain open.
No full configuration matrix or external-solver run was performed.
