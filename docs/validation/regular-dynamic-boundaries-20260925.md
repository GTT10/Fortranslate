# Regular Dynamic LES physical boundaries — 2026-09-25

Implementation `9fdd346`, based on `902a99f`.
`compute_dynamic_les_boundary_rhs_3d` connects regular-grid wide donors to the
existing ghosted Dynamic kernel. It supports reservoir inflow, extrapolated or
pressure-target outflow, slip and adiabatic/isothermal no-slip walls, including
prescribed tangential wall motion. The filter determines halo width (radius
plus five); repeated reflection supports domains smaller than that halo.
Physical-face flux integration returns the net inventory rate into the domain.
The padded-workspace maximum diffusivity remains a conservative bound.
RHS, diffusion maximum and optional boundary rate publish only on success.

`advance_dynamic_les_ssprk2_3d` appends optional physical inputs and a boundary
change output. Both stages enforce the existing diffusion bound. State,
temperature and integrated boundary change publish only after final recovery;
zero dt validates the first-stage operator and returns zero boundary change.
The periodic RHS path and existing weighted-RHS SSPRK2 arithmetic are retained.
The regular spray coupling/application still has its periodic-only Dynamic
guard; connecting this API there and extending the regular MPI path remain.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug. Fixed tests-enabled unit/runtime and selected
H2O2 tests-disabled spray runtime built. Selected execution and MPI execution
were not run. Final `spray_les_dynamic_boundaries`: **1/1 passed**, 4.47 seconds
(total 4.51), no skips.

Four physical combinations cover reservoir inflow/outflow, slip plus moving
adiabatic walls, moving isothermal walls, and a single pressure-target outlet.
SSPRK2 state agrees with the existing AMR single-root operator within eight
machine epsilons times each component magnitude; temperature agrees exactly.
Boundary budgets agree to 1e-12 relative maximum budget. Integrated inventory
change agrees with boundary change to 2e-14 times initial inventory scale.
Excessive dt preserves state, temperature and boundary change. Invalid wall
input preserves RHS, diffusion and rate outputs. The periodic wrapper RHS and
diffusion bound equal the original periodic operator exactly, with zero budget.

The initial new test incorrectly demanded bitwise state equality with AMR and
failed in 0.97 seconds. A diagnostic rerun measured 2.9103830456733704e-11
maximum state difference in the first inflow case and zero temperature difference.
AMR averages two Euler states while the existing regular integrator adds a
weighted RHS, giving different floating-point grouping. The new comparison was
corrected to the componentwise machine-epsilon bound above; no existing
acceptance tolerance or implementation arithmetic was changed to obtain parity.
The final test passed all four combinations. Mechanism provenance, Python
syntax of the contract checker, project contract and whitespace checks passed.
This is operator integration evidence, not spray application execution, MPI
physical-boundary support, external solver parity or physical validation.

## Commands and artifacts

Logs: `build/regular-dynamic-boundary-checks/` (fixed/selected builds,
initial `ctest.log`, final `ctest-final.log`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_dynamic_boundaries --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_dynamic_boundaries$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
