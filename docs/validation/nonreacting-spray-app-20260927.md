# Reaction-free selected spray application — 2026-09-27

Implementation: `e83c0d72e42190a4b969e6c521da0960935ba47e`.
The selected serial spray front end now accepts an explicitly nonreacting
bundle. CMake checks its boolean marker and zero reaction count. The shared
runtime keeps its usual nonempty-reaction requirement unless a caller explicitly
opts into nonreacting preparation. The spray application rejects chemistry ON
with an empty reaction array before creating output.

The [heptane fixture](../../cases/spray_3d/heptane_nonreacting.nml) uses pinned
NC7H16/O2/N2 material data, a periodic 4³ grid, four parcels, 400 K gas and
310 K liquid, two `1e-7` s steps. The upstream 273 K gas/272 K liquid and 2D
legacy boundary input are not reused: those temperatures are below the
imported NC7H16 NASA7 range. Source coefficients/ranges were not modified.
The fixture's Antoine evaluation interval 300–371.6 K is a computational
restriction, not new evidence for the fit's physical validity.

## Fresh build and completed checks

GNU Fortran 13.3.0 Debug, serial, spray ON, tests/MPI/SUNDIALS OFF. Fresh build
directory `build/heptane-spray-debug`; selected source is the hash-verified
bundle from [the import record](nonreacting-spray-import-20260927.md).

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/heptane-spray-debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DCMAKE_BUILD_TYPE=Debug -DPELEF_ENABLE_SPRAY=ON -DPELEF_ENABLE_TESTS=OFF -DPELEF_ENABLE_MPI=OFF -DPELEF_ENABLE_SUNDIALS=OFF -DPELEF_MECHANISM_BUNDLE="$PWD/build/heptane-3sp-reference/heptane_3sp.json" -DPELEF_MECHANISM_SOURCE="$PWD/build/heptane-3sp-reference/heptane_3sp.yaml"
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-spray-debug --target pelef_spray_3d_selected --parallel 2
python3 tests/application/check_nonreacting_spray_app.py --exe build/heptane-spray-debug/pelef_spray_3d_selected --source . --work build/heptane-spray-checks
```

Build succeeded. The final checker completed one evaporation run and two
expected rejections (chemistry ON; gas temperature 273 K), with zero unexpected
failures/skips. Rejections produced no gas CSV. Results:

- Evaporated liquid mass: `4.516827337669613e-15` kg.
- Integrated gas NC7H16 mass: `4.516827337676076e-15` kg.
- Maximum normalized total-mass error: `1.0568263930649006e-15`.
- Maximum normalized total-energy error: `1.2270212629072912e-16`.

Gates: relative vapor/liquid transfer mismatch below `1e-9`; mass/energy
inventory errors below `1e-12`. No tolerance changes were made. Mechanism
provenance, Python syntax, project contract and diff checks passed.

Artifacts: `build/heptane-spray-checks/nonreacting_qxbb_f66/` (inputs, logs,
CSV/checkpoint output and summary). Initial direct run:
`build/heptane-spray-run/`. Logs: `build/heptane-spray-configure.log`,
`build/heptane-spray-build.log`. The initial configure exposed the old blanket
zero-reaction rejection; configuration succeeded after the explicit contract
was implemented. No full CTest matrix or full PeleC spray comparison was run.

This establishes reaction-free application integration and bounded exchange
conservation. It is not a reproduction of the original 273 K evaporation case
or experimental evaporation validation. Other selected front ends retain
their existing reaction-count requirements.
