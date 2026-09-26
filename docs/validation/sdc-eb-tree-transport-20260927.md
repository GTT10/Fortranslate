# Provisional EB tree SDC transport — 2026-09-27

## Source and behavior

Implementation `3c9d8fdb`, based on `096846f6` on local `main`.
PeleC pin remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

`compute_sdc_eb_tree_transport_3d` evaluates the existing molecular transport
stage on every EB patch. It returns packed rates per fluid volume,
time-integrated face flux densities without aperture, and the minimum limiter
coefficient. State and temperature remain inputs. Gas halos retain provisional
parent values through `restrict_snapshot=.false.`; no averaging, redistribution
or reflux occurs here. Rates, fluxes and limiter publish together after all
patches succeed. Existing physical-face and wall-temperature dispatch is reused.

## Completed focused verification

GNU Fortran 13.3.0, Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`. Existing compiler warnings remain.

- Mechanism provenance preflight passed for `mechanisms/h2o2_full.json`.
- `test_eb_tree_transport_3d` built successfully.
- Final CTest `unit_eb_tree_transport_3d`: 1/1 passed, 4.38 s test / 4.43 s total.
- New SDC check: z-plane cut, three levels, five patches, viscosity,
  conduction, diffusion and barodiffusion enabled. Reconstructing old state
  plus dt times rate and applying existing flux synchronization reproduces
  the common Euler transport state within 2e-12 times max(1, abs(reference)).
  The limiter agrees exactly.
- A negative fluid-cell density rejects without changing any packed rate,
  x/y/z flux array or limiter. Existing transport boundary and conservation
  checks remain in the same executable.
- The initial rejection test failed because a negative temperature guess was
  recoverable from the conserved state. It was changed to an invalid density;
  no numerical tolerance or solver acceptance rule was relaxed.
- Python syntax, project contract and `git diff --check` passed.

Commands from repository root:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_transport_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_transport_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Logs: `build/sdc-eb-tree-transport-checks/`, including the initial failed check.
No coupled EB tree SDC trajectory, selected build, MPI run, full matrix,
installed-binary audit or physical validation was performed. Selected-level
transport and ancestor time interpolation are not connected by this interface.
