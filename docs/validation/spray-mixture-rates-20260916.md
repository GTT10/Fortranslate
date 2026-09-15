# Mixed-liquid rates and frozen-gas transport — 2026-09-16

Tested parent: `eabf686` plus this change. Exact modified source/test hashes and
logs are retained in `build/mixture-rates-checks/`.

## Implementation and source boundary

The film builder uses Raoult surface mole fractions, suppresses components with
no outward vapor driving force, and forms a one-third-rule film state. Component
binary diffusivities are evaluated with the existing harmonic carrier-mole average.
Their surface-mole-weighted mean supplies the common transfer correlation;
component mass rates divide in proportion to X_surface*D. This weighting structure
is recorded in the pinned Drag.H source in `references/spray_mixture.json`.

Drag and selectable Ranz–Marshall/Abramzon–Sirignano heat/mass transfer reuse the
existing correlations at mixed-liquid density and cp. The temporary property
record clears unused fit coefficients after evaluation. Component enthalpies and
energy remain resolved individually. A frozen-gas interval now combines drag,
flight, heating, evaporation and returned conservative gas/species increments.
It rejects overdraw or invalid thermal states without changing parcel/component
storage. Zero duration preserves the input exactly.

PeleF's carrier molecular-weight, direct binary-diffusion and film-cp heat
conventions are retained. No exact upstream coefficient or timestep parity is
claimed. The caller still owns adaptive interval selection and grid feedback.

## Focused verification

GNU Fortran 13.3.0 Debug, MPI disabled, bounds/FPE checks:

```sh
cmake --build build/port-spray-debug --target test_spray_core --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(mixture_rates|mixture_exchange|mixture)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

All three tests passed in 0.03 s with no failures or skips. Provenance, contract
and whitespace checks passed; no Python source changed. The first constant-property
run passed; review then identified stale polynomial coefficients in the temporary
constant property record. Clearing them and extending the case to polynomial
density/cp also passed the final focused run.

The rate test covers both transfer-correlation choices in the single-component
limit: drag, heat and evaporation agree with the existing pure-liquid route to
2e-13 relative tolerance. The two-component case checks preferential evaporation,
polynomial liquid properties, drag/flight, total mass/momentum/energy conservation,
gas species closure, zero-duration identity, full rollback for an oversized
interval, and suppression of one component by its ambient vapor. H2O/H2O2 values
are illustrative fixtures, not physical mixture validation data.

## Remaining work

The mixed-liquid path is still an API using an external component-mass array.
Persistent particle composition, regular-grid feedback, mixture wall/breakup
handling, CLI selection and coupled restart remain open. No full configuration
matrix, external trajectory comparison or physical validation was run.
