# Two-level EB hydro viscosity — 2026-09-25

Code `72f47ce`, parent `6051259`. Common and subcycled two-level Euler hydro
accept optional finite nonnegative `difmag`, default zero. Positive values
reach coarse and fine Euler fluxes before averaging/reflux. Fine recursive
steps rebuild complete halos from time-interpolated parent conserved states.
Transport-only use of the private dispatcher rejects positive coefficients.

The patch halo API optionally returns volume fractions and complete interior,
face, edge and corner states. Physical-domain indices are clamped; interior
uses fine data and exterior patch donors use PCM parent data. Previously checked
open-face values are retained. Covered storage stays zero without EOS recovery.
State, temperature and volume fractions publish together only after success.
Omitting the optional argument retains the existing face-only path.

PCM ancestor volume fractions remain an approximation, not reconstructed
fine ghost geometry. Multipatch peer halo overlay and full reaction/transport
wrapper coefficient inputs are not connected by this change. Application
input wiring and the uniform reference default policy remain open here.

## Completed checks

GNU Fortran 13.3.0, Python 3.12.3, Ninja Debug; reused
`build/spray-amr-reactive-debug`. CTests **2/2 passed**, default 0.97 s and
coefficient 0.037 1.13 s, total 2.14 s. No failures or skips. Compiler warnings
remain. Provenance, Python syntax, contract and whitespace checks passed.

Ratios 2/3, interior/physical-boundary patches, common/recursive stepping and
three successive steps satisfy existing 2e-12 composite inventory/reflux bounds.
Existing timestep/geometry/ratio rollback tests pass. Positive mode also checks
complete donor states and volume fractions at every halo index (including
edges/corners), using floor-based parent indexing as the expected mapping.
A bad fine donor rejects with all three halo outputs unchanged. These tests do
not establish full PeleC geometry parity or physical validation.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_two_level_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_two_level(_viscosity)?_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/eb-two-level-viscosity-checks/` contains `build.log`,
`rebuild.log` and `LastTest.log`.
