# Nozzle Mach inlet state construction — 2026-09-25

Implementation `fcb9444`, based on `1294e04`.
`mach_inlet_ghost_3d` ports the pinned EB-ConvergingNozzle inlet-state equations,
generalized from its low-x face to all six Cartesian orientations. It preserves
the interior mass fractions, uses interior pressure by default or an explicit
pressure target, and recomputes density/internal energy at prescribed temperature.
The inward normal speed is Mach times the sound speed of that inlet state.
Optional supplied velocity fluctuations are added in Cartesian coordinates;
optional total velocity replaces the mean. The two velocity options are mutually
exclusive. These are supplied vectors, not a turbulent-inflow generator.

The API accepts Mach 0..100, positive EOS-valid temperature and positive pressure
up to 1e15 Pa. It publishes ghost state and temperature only after successful
EOS reconstruction. Invalid face, shapes, finite/range inputs, thermodynamics
or conflicting velocity modes reject. It does not classify the flow Mach regime,
provide GC-NSCBC, or yet connect the inlet mode to solver ghost filling/input.

## Frozen reference

PeleC commit `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
Local reference artifacts were rehashed against the identities recorded in
[the pressure-ghost source record](characteristic-pressure-ghost-20260918.md):

- `build/nscbc-reference/nozzle-prob.H`, SHA256
  `8e3d840a7857ec4d765019aad00369ea6fbe9aef5731fc3970c79c29e97c0bcb`.
  Lines 84–136 implement interior composition, inlet types 0/1 for pressure,
  prescribed temperature/Mach, optional fluctuating/total velocity, and conserved
  state construction. The comment mentions inlet type 2 but the code rejects it;
  this implementation does not invent a type-2 equation.
- `build/nscbc-reference/BoundaryConditions.rst`, SHA256
  `897b100e6ab931a850d0427dc324852ef3dfc98879772ba0c4c79351609c4f2d`.
  The guide explicitly says GC-NSCBC was not ported from Fortran into C++.
  This kernel therefore ports the implemented nozzle inlet, not deprecated NSCBC.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug. Fixed core/unit and selected H2O2 spray runtime
built. Selected execution was not run. `unit_reactive_boundaries_3d`:
**1/1 passed**, 5.89 seconds, total 5.93. No failures or skips.

Six faces and three modes verify inward Mach velocity, supplied fluctuations,
total-velocity replacement, temperature, pressure, ideal-gas density ratio and
frozen composition at a 3e-12 scaled bound. Negative pressure/temperature and
conflicting velocity options preserve sentinel outputs. Existing physical hydro,
pressure-outlet and reservoir-boundary checks remain in the same test.
No existing tolerance was changed. Mechanism provenance, Python syntax of the
contract checker, project contract and whitespace passed. This is state-kernel
evidence; no full C++ case comparison, solver integration or physical validation
was run.

## Commands and artifacts

Logs: `build/mach-inlet-checks/`.

```sh
sha256sum build/nscbc-reference/nozzle-prob.H build/nscbc-reference/BoundaryConditions.rst
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_boundaries_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_boundaries_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
