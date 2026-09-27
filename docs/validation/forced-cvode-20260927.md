# Forced CVODE reactor — 2026-09-27

Implementation `ac0ecd1d`, based on `0d0d600e`, local main, version 0.246.0.

## Change

`sundials_forced_reactor_mod::advance_forced_reactor_cvode` integrates all
species densities with CVODE BDF, using `forced_reactor_rhs` for chemistry,
constant external species rates and a linearly changing internal-energy
density. CVODE forms a dense numerical Jacobian and enforces nonnegative
species densities. A stop time prevents evaluation beyond the source interval.
The complete interval publishes density, temperature and step count only after
the solve, final thermodynamic recovery and resource cleanup succeed.

Private synchronous callback slots isolate the forcing data from the existing
constant-volume CVODE registry. Callback tokens remain alive until CVODE memory
is freed; failed solves release their slot. Concurrent threaded calls remain
unsupported, consistent with the existing backend. The module is included in
selected and spray CVODE runtime source lists.

## Independently completed checks

Fresh Ninja Debug build `build/forced-cvode-debug`, GNU Fortran 13.3.0 and local
SUNDIALS 7.2.0 official Fortran modules; MPI and spray disabled for this focused
backend check. Configure uses the pinned H2O2 source/bundle. Artifacts:
`build/forced-cvode-checks/` (configure/build/provenance logs and initial/final
CTest/LastTest logs).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/forced-cvode-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DPELEF_ENABLE_MPI=OFF -DPELEF_ENABLE_SUNDIALS=ON -DPELEF_ENABLE_TESTS=ON -DPELEF_MECHANISM_BUNDLE=/home/gtt13/projects/pelec/mechanisms/h2o2_full.json -DPELEF_MECHANISM_SOURCE=/home/gtt13/projects/pelec/mechanisms/h2o2_cantera.yaml -DSUNDIALS_DIR=/home/gtt13/projects/pelec/build/sundials-7.2-fortran/install/lib/cmake/sundials
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/forced-cvode-debug --target test_sundials_forced_reactor test_sundials_constant_volume_reactor test_sundials_constant_volume_multi_context --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/forced-cvode-debug -R '^unit_sundials_(forced_reactor|constant_volume_reactor|constant_volume_multi_context)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/forced-cvode-debug -R '^unit_sundials_forced_reactor$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The final focused test passed in 0.13 s (CTest total 0.16 s). Existing
constant-volume and multi-context tests passed in the first invocation,
0.01 s each. That invocation failed the new test because the native comparison
solver could not complete the 1e-7 s interval at the requested 1e-10/1e-13
tolerances. The final test retains CVODE mass/energy checks over 1e-7 s and
compares both solvers over 1e-10 s with those same tolerances. The relative
species difference is 3.7993516232764078e-11, below the unchanged 2e-7 criterion;
temperature also meets the unchanged 2e-4 K criterion. Longer-interval
native/CVODE agreement is not established.

Seventy inert forced solves agree with the linear species-density solution and
energy-source balance, exercising slot reuse beyond registry capacity. An
intentional one-step limit rejects a reacting solve without changing caller
state and allows a subsequent solve. The expected CVODE `mxstep` diagnostic
appears in the passing test log. Build, provenance, Python syntax, contract and
whitespace checks passed; existing real-comparison warnings remain.
Pinned mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Remaining integration

`react_sdc_state_3d` still invokes the native forced reactor. Backend dispatch
and propagation of the selected integrator through SDC applications are required
next. This library increment does not establish SDC/CVODE application behavior,
MPI or spray builds, detailed-fuel performance, installed-binary behavior,
external full-solver parity or physical validation.
