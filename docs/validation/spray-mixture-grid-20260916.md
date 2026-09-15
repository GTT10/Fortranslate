# Mixed-liquid grid feedback and composition restart — 2026-09-16

Tested parent: `c196ed8` plus this change. Exact modified source/test hashes and
logs are retained in `build/mixture-grid-checks/`.

## Implementation

The regular-grid mixture API retains component masses per droplet in an allocated
component-by-parcel matrix. It samples gas with the existing stencil, bounds
component/thermal/drag/CFL changes and retries inadmissible feedback with a shorter
interval. Only a complete accepted update replaces gas, temperatures, parcels,
composition and boundary budgets. Component inventory is included in system totals.

Periodic, reflected and open boundaries reuse the existing geometric handling
with mixture energy/density. Open boundaries remove all component masses and retain
their species budgets. ETAB/KHRT representative breakup uses prescribed mixture
viscosity/surface tension and mixed density; component masses scale with droplet
mass. Computational splitting copies per-droplet composition into descendants.

Optional checkpoint composition storage adds `_MIX_V1`, a bounded component count
and per-parcel component-mass records. Strict finite/nonnegative/mass-closure
validation and atomic reader commit preserve the previous state on malformed input.
The caller's context binds the ordered component property table and settings.

## Focused verification

GNU Fortran 13.3.0 Debug, MPI disabled, bounds/FPE checks enabled:

```sh
cmake --build build/port-spray-debug --target test_spray_core test_spray_integration --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_integration_(mixture_grid|checkpoint)$' --output-on-failure
ctest --test-dir build/port-spray-debug -R '^spray_integration_mixture_grid$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

The two-entry run passed in 0.09 s. After adding malformed-composition read and
full particle-state restart assertions, the mixture entry passed in 0.06 s.
No failures or skips occurred. Provenance, contract and whitespace checks passed;
no Python source changed.

The 5³-cell case checks evaporation into two gas species; component, total mass,
momentum and energy conservation; continuous-versus-restarted equality of gas,
temperature, parcel position/velocity/temperature/mass and component arrays;
rejection when composition output is omitted; malformed component-row rollback;
late substep-limit rollback; mixed outflow species/energy budgets; and representative
KHRT breakup plus computational splitting with retained component inventory.
The old pure checkpoint test also passed. H2O/H2O2 properties are synthetic fixtures,
not physical mixture qualification. Reflection is implemented through the shared
boundary routine but was not separately exercised in this focused increment.

## Remaining work

This is a grid-coupling/checkpoint API; mixture namelist input, scheduled injection
and the combined hydro/chemistry application driver are not yet connected. Mixture
wall impact/films and KH shed-mass children remain open. No particle MPI/AMR/EB,
full matrix, external trajectory parity or physical validation is claimed.
