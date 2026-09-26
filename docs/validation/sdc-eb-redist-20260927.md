# Small-cell redistribution as an SDC source rate — 2026-09-27

Implementation `0eace65f`, base `2c43bba2`; PeleF 0.246.0.
PeleC pin `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

## Implementation

`redistribute_sdc_eb_rate_3d` takes geometry, old state, conservative rate,
positive dt and `none`, `flux_redist` or `state_redist`. Flux redistribution
acts directly on rates using the existing neighbor algorithm. State
redistribution acts on old state plus dt times rate using the existing weighted
state algorithm, then returns the difference from old state divided by dt.
Covered inputs are excluded and covered output rates are zero. Invalid active
inputs or failed redistribution leave the public rate unchanged. No temperature
recovery or chemistry is performed here; generic signed components are supported.

Pinned `Source/Hydro.cpp` lines 720–798 were inspected with `git show`: divergence
is assembled before `ApplyMLRedistribution`, whose output feeds postprocessing.
This motivates the source-rate interface. The existing local algorithms are not
claimed identical to all AMReX redistribution modes. Reference density-from-
species postprocessing and optional mass-fraction cleanup are separate from this
helper. It is not yet connected to a full EB SDC timestep.

## Focused evidence

Artifacts: `build/sdc-eb-redist-checks/` (local). GNU Fortran 13.3.0,
existing Ninja Debug tests-enabled build `build/spray-amr-reactive-debug`;
compiler/backend unchanged.

- Provenance passed: `python3 tools/check_mechanism_provenance.py --bundle
  mechanisms/h2o2_full.json`; SHA
  `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
- Build passed: `env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build
  build/spray-amr-reactive-debug --target test_eb_reactive_redistribution_3d
  --parallel 2` (`build.log`). Compiler warnings remain.
- CTest **1/1 passed**, 0.01 s (0.06 s total): `env -u LD_LIBRARY_PATH
  -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R
  '^unit_eb_reactive_redistribution_3d$' --output-on-failure` (`ctest.log`).
- Added an 8x1x1 plane-cut fixture with covered NaNs, nonuniform old fields and
  signed rates. All three modes conserve volume-weighted rate integrals within
  2e-12 of the test scale. The state-mode rate reconstructs the independently
  invoked weighted redistribution result within 1e-13. Covered rates are zero.
  Invalid active rate rejects while retaining the previous output exactly.
- Existing axis, transverse-neighbor, reactive-update, curved and connected
  neighborhood checks also passed. Python syntax for
  `tools/check_project_contract.py`, project contract and `git diff --check` passed.

No coupled EB SDC, application trajectory, MPI path, selected build, full matrix,
installed audit or physical validation was run. State-to-rate subtraction has
ordinary floating-point cancellation at very small dt; no broad timestep
accuracy claim is made. The full port remains incomplete.
