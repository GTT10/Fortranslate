# Regular spray Dynamic LES physical boundaries — 2026-09-25

Implementation `0ab5a3e`, based on `886386e`.
The serial coupled spray API and fixed/selected application now accept the
existing physical gas/particle boundaries with Dynamic LES. The two Dynamic
half steps receive reservoir, wall and pressure inputs, and return exchange
into the existing coupled boundary ledger. Candidate gas, temperature, parcels,
mixture masses and films retain the existing whole-step publication contract.
The periodic-only guards were removed after the regular Dynamic operator gained
physical donors and stage budgets. Positive cs enables Dynamic SGS; zero still
disables it. Filter settings and physical inputs retain checkpoint identity.
Existing CFL selection and stage-failure timestep halving remain unchanged.

The boundary checker adds opt-in `--dynamic-physical`, preserving its default
periodic Dynamic fixture. Changing a periodic boundary at restart is now an
incompatible-checkpoint error rather than an unsupported-model error. Existing
physical cases already check wall/inlet changes. Regular MPI Dynamic physical
boundaries remain unconnected; this change concerns the serial regular path.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug. Fixed tests-enabled and selected H2O2
tests-disabled spray applications built. Selected execution was not run.
CTest `spray_dynamic_wall_restart`: **1/1 passed**, 37.67 seconds (total 37.71).
No failures, skips or numerical tolerance changes.

The fixture uses Gaussian test filtering, test ratio 4, grid ratio 2, a 4-cubed
domain, heated isothermal walls and one prescribed tangential wall speed.
Two parcels remain live. Continuous and checkpoint-resumed gas/parcel CSVs
are byte-identical; final history rows, including boundary ledgers, are equal.
Net wall heat is positive. Existing mass, momentum and energy residual checks
pass at their 1e-9 scaled bounds. Changed wall motion, model/filter controls
and invalid settings reject without final histories. Inflow, pressure outlets,
other closures and selected execution were not rerun in this application check.
The previous operator test supplies their narrower operator evidence.

Mechanism provenance, Python syntax for the changed checker and project-contract
checker, project contract and whitespace passed. This is focused coupled/restart
evidence, not MPI physical-boundary execution, physical validation, external
solver parity or a full configuration matrix.

## Commands and artifacts

Logs: `build/regular-dynamic-spray-checks/`.
Application artifacts:
`build/spray-amr-reactive-debug/tests/spray_dynamic_wall_restart/spray-boundaries-9mshorto/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_spray_3d_selected --parallel 2
ctest --test-dir build/spray-amr-reactive-debug -R '^spray_dynamic_wall_restart$' --output-on-failure
python3 -m py_compile tools/check_spray_boundaries.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
