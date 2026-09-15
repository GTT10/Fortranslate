# KH shedding with computational splitting — 2026-09-16

Tested parent: `f6e47cf` plus this change. Exact modified source/test hashes
and logs are retained in `build/khrt-split-checks/`.

## Implementation

KH shedding can now be combined with a positive `maximum_multiplicity`.
Computational splitting divides both initial drop count and accumulated N*r³
among the pieces, rather than copying extensive state into every descendant.
The final piece carries the floating-point remainder. Mass, momentum, energy,
lineage and existing atomic allocation/rollback rules retain their meanings.
No checkpoint format change is needed: the existing KH state records contain
these partitioned values. This is a PeleF computational representation choice;
no additional upstream physical model was introduced.

## Focused verification

GNU Fortran 13.3.0 Debug, MPI disabled, bounds/FPE checks enabled:

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d test_spray_core --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(split|khrt)$' --output-on-failure
python3 tests/spray/test_spray_application.py --exe build/port-spray-debug/pelef_spray_3d SprayApplication.test_khrt_children_and_accumulation_restart
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Both CTest entries and the selected Python test passed without failures or
skips. The application test covers four configurations: release and
accumulation-only, each with and without computational splitting. It checks
multiplicity caps, system mass/momentum/energy conservation, retained nonzero
accumulators and exact restart CSV/checkpoint bytes. The split kernel test
checks repeated division with remainders, conserved total count and shed volume,
per-drop ratios, and preservation of KH state on rejected capacity limits.
Provenance, syntax, contract and whitespace checks passed.

## Boundary

The global jet normalization and per-parcel release threshold remain unchanged;
changing the computational multiplicity can change release timing. Conservation
and restart equality for a fixed setting do not prove trajectory invariance
across partition choices. KH children with wall impact, general turbulent
dispersion, mixture droplets and particle MPI/AMR/EB remain open. No full test
matrix or physical validation was run.
