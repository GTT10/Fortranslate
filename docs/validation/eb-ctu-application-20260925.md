# Single-level EB CTU application integration — 2026-09-25

## Implementation and scope

Code `350ad07`, parent `4feb904`, adds `hydro_method='plm_ctu'` to the fixed
and configure-time selected single-level reactive EB applications. The default
is `pcm`. CTU requires `riemann_solver='pelec'`; other applications using the
shared reader reject this selector unless explicitly enabled.

The application builds six-cell extended geometry once. Both chemistry/hydro
Strang and full chemistry/transport/hydro composition accept this optional
geometry and dispatch the outflow CTU entry. The central geometry is checked
against the physical domain (bounds, classification, all cell/face centroids,
fractions and wall metrics) before chemistry runs. CTU retains fourth-order
slope selection, multidirectional flattening and the existing viscosity and
redistribution controls. The outflow wrapper restores original covered-cell
storage after an accepted update, satisfying the pre-existing split-operator
contract while keeping covered internal hydro storage inactive.

Checkpoint solver identity appends `hydro-v1 plm_ctu` for the new route;
PCM identities retain their previous form. The minimal numerical example is
`cases/reactive_eb_3d/plm_ctu.nml`. Reference algorithm pin remains
PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

## Completed checks

GNU Fortran 13.3.0 / OpenMPI 4.1.6 / Python 3.12.3. Existing Ninja Debug builds:
`build/spray-amr-reactive-debug` (tests enabled) and
`build/spray-amr-reactive-selected-debug` (tests disabled, selected
`mechanisms/h2o2_full.json`). No install-prefix or compiler/backend changes.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_eb_3d test_reactive_eb_3d_config --parallel 2
python3 tools/check_eb_ctu_application.py --executable build/spray-amr-reactive-debug/pelef_reactive_eb_3d --work-dir build/eb-ctu-application-checks/app
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected --parallel 2
python3 tools/check_eb_ctu_application.py --selected --executable build/spray-amr-reactive-selected-debug/pelef_reactive_eb_3d_selected --work-dir build/eb-ctu-application-checks/selected-app
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_eb_3d_config$' --output-on-failure
python3 -m py_compile tools/check_eb_ctu_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Both application checker invocations passed. Each checks five successful runs
(inert hydro, reacting hydro, coupled reference, checkpoint stop, resume), plus
three expected rejections (method mismatch at restart, incompatible solver,
unknown method). Coupled restart CSV matches uninterrupted output byte for
byte; chemistry and transport each change the result. CSV values are finite,
species sum to unity, and per-component inventory ledgers satisfy the checker.
The existing application inventory check also runs internally. These are two
focused application scenarios, not a configuration matrix or external parity.
The existing configuration unit passed **1/1**, 0.06 s / 0.10 s total.
The checker is registered as `reactive_eb_3d_ctu_application` for future CTest
runs; the two executions above used its direct Python entry.

The initial fixed inert trial failed the existing covered-cell preservation
check. The wrapper was corrected to restore those values and both complete
checker invocations then passed. No tolerances or golden hashes changed.
There were no build failures or skips. Compiler warnings remain.

Successful artifacts:
- `build/eb-ctu-application-checks/app/run-on9_6cp6/`
- `build/eb-ctu-application-checks/selected-app/run-1jn1bgnk/`
- build logs and `LastTest-config.log` in `build/eb-ctu-application-checks/`.

The initial failed trial is also retained under that directory's `app/`.
Provenance matched SHA
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Python syntax, project contract and diff checks passed.

## Remaining boundary

The runtime scenarios use plane geometry and H2O2. Extruded/nodal CTU application
runs, other physical boundary donors, EB tree/MPI/spray composition, broad
convergence/scaling and physical validation remain open. No installed-binary
audit, external solver run or full build matrix was performed.
