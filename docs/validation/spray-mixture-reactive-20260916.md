# Mixed-liquid reacting-flow integration — 2026-09-16

Tested parent: `a91eb3f` plus this change. Exact modified source/test hashes and
logs are retained in `build/mixture-reactive-checks/`.

## Implementation

The existing symmetric spray/LES/chemistry/hydro driver accepts the optional
mixture property table and component-mass matrix together. Its two spray stages
use the mixture grid update; hydro, molecular transport, LES and chemistry retain
their existing implementations and ordering. External acceleration continues to
account for gas/particle impulse and work.

Composition is copied with the candidate gas and parcel state. A later failed
stage discards all candidate state and returns zero boundary/source ledgers.
Only complete acceptance commits the mixture mass matrix. Missing property/storage
pairs, simultaneous per-port pure-property tables, particle wall impact and
nonempty films reject explicitly. Existing pure-liquid call signatures retain
all preceding positional argument meanings.

## Focused verification

GNU Fortran 13.3.0 Debug, MPI disabled, bounds/FPE checks enabled:

```sh
cmake --build build/port-spray-debug --target test_spray_integration --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_integration_(mixture_reactive|coupled|rollback)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

All three tests passed in 24.72 s: existing pure coupled 4.69 s, pure rollback
0.10 s and mixed reactive 19.91 s. No failures or skips. Provenance, contract and
whitespace checks passed; no Python source changed.

The mixture test uses 5³ cells, reacting H2/O2/N2 gas initially at 1400 K,
illustrative H2O/H2O2 droplets, native implicit chemistry, molecular transport,
Smagorinsky coefficient 0.1, and acceleration -9.81 m/s² in z. Both liquid
components decrease and gas H2 changes. Mass and total momentum/energy close
with the external-force ledger. Stop/restart yields exact gas, temperature,
component mass, parcel position/velocity/temperature/mass and force-output equality.
An invalid chemistry backend encountered after the first mixture stage restores
all these input states and zeroes returned force. The missing property-table
case also rejects. These synthetic properties are implementation fixtures.

## Remaining work

Mixture CLI inputs and scheduled injection are not connected yet. Mixture wall
impact/films, KH shed children, particle MPI/AMR/EB and physical validation remain
open. The new focused run does not qualify CVODE or a full configuration matrix,
and does not establish external solver or physical ignition parity.
