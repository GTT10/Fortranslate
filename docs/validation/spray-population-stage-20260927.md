# Combined EB level population stage — 2026-09-27

Implementation `c962220a`, based on `ad3eab58`, local `main`, version 0.246.0.
Pinned PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Implementation

`advance_mpi_eb_population_stage` accepts three level populations ordered
physical, virtual and ghost. It fills the stage-time EB gas halo once, with
optional parent state/temperature endpoints and interpolation fractions, then
advances all populations privately against that same frozen gas.

Physical and virtual populations use retained level routing; ghost populations
use grown-box filtering. Root ghosts must be empty. Each population contributes
same-level time-integrated gas feedback and an unmatched halo remainder. The
three source fields are added with collective overflow rejection. Only physical
wall/outflow ledgers are exported as physical-system inventories. Gas inputs
remain read-only. All populations, combined sources and ledgers publish only
after all three updates and source sums succeed.

Transfers occur per population followed by field addition. This is a linear
source assembly with a different floating-point addition order from one shared
temporary source and one transfer; no bitwise reference parity is claimed.
The caller remains responsible for creating compatible populations, their
old/new SDC lifecycle, breakup/children/films and external particle acceleration.
Complete recursive application dispatch remains open.

## Completed checks

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6. `spray_mpi_eb_sample_2` and `_4` passed:
2/2, 0.64/0.63 s; CTest total 1.31 s. No build/test failures or skips occurred.
Provenance, Python syntax, project contract and whitespace checks passed.

The nonperiodic planar EB fixture places one particle from each population at
the same fine ghost-cell center with distinct IDs. Parent fractions 0, 0.25 and
1 give all three exactly the same velocity, temperature and component updates
as the separately checked ghost stage. Their combined remainder equals three
times the single-population remainder within 1e-12 relative tolerance; all
valid-cell sources are zero in this ghost-center fixture. The independent
frozen-gas kernel comparison and fluid-volume source inventory checks still run.

An invalid node in the third population rejects after the first two private
updates. Component masses, nodes, velocities, temperatures, combined sources,
remainders and sentinel ledgers remain unchanged. This explicitly exercises
late failure, not only an invalid initial gas halo.

The new combined fixture has zero drift and no wall/outflow event; nonzero
combined valid-cell deposition, root-level execution and copy-boundary ledger
exclusion are implemented but not independently exercised here. This is not
complete recursive population evolution, a full configuration matrix, an
installed audit, external solver numerical parity or physical validation.

## Reproduction

Artifacts: `build/spray-population-stage-checks/`, including build, provenance,
contract, CTest summary and copied detailed `LastTest.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_sample --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sample_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
