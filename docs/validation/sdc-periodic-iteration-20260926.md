# Nonreacting periodic SDC iteration — 2026-09-26

Implementation `ee4a46fd`, parent `900a8d25`, local main, version 0.246.0.

`advance_reactive_sdc_periodic_3d` evaluates old source rates once and initializes
new rates from them. Each iteration forms the half-new/half-old predictor,
recomputes CTU fluxes from the original state, assembles the candidate, recovers
temperature, evaluates new rates at time+dt, and rebuilds the state from the
original state and updated source list. The callback receives physical cells;
periodic source halos are supplied by CTU. This follows the nonreacting branch
of pinned `Advance.cpp::do_sdc_iteration` and `Hydro.cpp` source preparation.
Any source, flux or EOS failure preserves the caller's entire state/temperature.

## Verification

Final `unit_reactive_ppm_step_3d` passed **1/1**, 5.27 s test / 5.32 s total.
New uniform fixtures use two state-dependent source families, +12U and -2U.
One and three iterations agree with their analytic trapezoidal fixed-point
recurrences, and source callbacks are counted (one old evaluation plus one
new evaluation per iteration). A failed new-source callback verifies rollback.
Existing CTU/source/split checks also ran. No test failures or skips.

An initial successful run took 5.39 s test / 5.44 s total. Its internal test
callback caused an executable-stack linker warning. Moving that callback to
a test module removed the warning; the final rerun above uses that version.
Provenance, Python syntax, project contract and whitespace passed. GNU Fortran
13.3.0, Python 3.12.3, existing tests-enabled Ninja Debug build; existing compiler
warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_step_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_step_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/sdc-periodic-iteration-checks/` contains provenance, initial
and final build/CTest logs and `unit.log`. Reference inspected via `git show`
at `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8` in
`/mnt/sn850/Projects/gtt2/external/PeleC`.

This is a library driver with focused analytic checks. Forced chemistry,
physical-boundary/AMR/MPI iteration, application selection, external executable
parity and physical validation remain open.
