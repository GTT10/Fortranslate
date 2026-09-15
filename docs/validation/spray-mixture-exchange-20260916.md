# Mixed-liquid surface and component exchange — 2026-09-16

Tested parent: `083dff4` plus this change. Exact modified source/test hashes and
logs are retained in `build/mixture-exchange-checks/`.

## Implementation

The surface API uses the pinned Drag.H liquid mole-fraction conversion and
Raoult relation. The existing saturation routine now accepts an optional liquid
mole fraction, applied after evaluating the actual-pressure surface cap. The
pure-liquid default remains one. A dilute component can therefore contribute a
valid partial pressure even when its pure saturation ratio is one or greater.
The combined surface must remain below saturation. Source identity is recorded
in `references/spray_mixture.json` from the pinned local Git blob.

The exchange API updates component masses per physical droplet from supplied
nonnegative evaporation amounts. Remaining composition determines liquid energy
and recovered temperature. Returned gas inventory maps each component to its
vapor species and includes evaporated momentum and total energy. Heat input is
for the whole parcel. Complete evaporation returns the entire finite inventory
without an undefined zero-mass energy inversion. Failures preserve both parcel
and component array and return zero gas exchange.

## Focused verification

GNU Fortran 13.3.0 Debug, MPI disabled, bounds/FPE checks enabled:

```sh
cmake --build build/port-spray-debug --target test_spray_core --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(mixture|mixture_exchange|surface_cap|antoine)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Four tests passed in 0.04 s with no failures or skips. Provenance, project
contract and whitespace checks passed; no Python source changed.
Tests verify mole-based surface partitioning, diluted pure saturation, total
saturation rejection, separate component losses and species routing, mass,
evaporated momentum and energy balances, component depletion, terminal evaporation,
overdraw rollback and the existing pure-property regressions. The preceding
mixture test still checks nonlinear cp and energy inversion. Illustrative H2O/H2O2
property values are synthetic fixtures, not physical validation data.

## Remaining work

Evaporation amounts and vapor enthalpies are supplied to this API. Mixture film
transport/rate calculation, ambient vapor suppression, drag/flight coupling,
persistent component storage, CLI selection and coupled restart remain open.
Many-to-one liquid-to-gas species mapping is rejected. This does not establish a
complete mixed-droplet CFD solver, full matrix qualification or physical validation.
