# SDC prescribed outlet pressure — 2026-09-26

Implementation `05688620`, parent `d248a9c7`, local main, version 0.246.0.

SDC now forwards optional outlet pressure through existing pressure ghost
filling, halo CTU boundary Riemann flags and physical transport fluxes. The
shared application driver passes its configured pressure array. Library
validation rejects nonfinite/negative/excessive pressure and positive pressure
on non-outflow faces. Input no longer rejects this supported SDC combination.

## Completed checks

- GNU Fortran 13.3.0 Debug fixed tests-enabled build
  `build/spray-amr-reactive-debug`, target `pelef_reactive_3d`, passed.
- **2 final application cases passed**: a 4³ nonuniform reacting hotspot with
  transport and external force, one step to 1e-9 s, two SDC iterations.
  Compared zero-gradient outflow with 290 kPa prescribed on x+ (initial
  interior pressure 300 kPa). Maximum x-velocity difference is
  **0.009872140578824684 m/s**. The pressure case's normalized Euler inventory
  residual is **1.25624e-17**, transport theta 1.
- Initial 250-kPa trial was rejected at pressure ghost filling. Existing
  extrapolation uses p_ghost = p_interior + 2*layer*(p_outlet-p_interior),
  producing negative pressure in the fourth halo layer for that input.
  Diagnostic instrumentation located the failure and was removed before the
  final build/runs. The final fixture uses admissible 290 kPa; no tolerances or
  ghost-state rejection rules changed. Rejected input/log and debug logs remain.
- Provenance, Python syntax, project contract and whitespace passed. No final
  failures or skips; existing compiler warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_3d --parallel 2
python3 build/sdc-pressure-outlet-checks/check.py
python3 -m py_compile tools/check_project_contract.py build/sdc-pressure-outlet-checks/check.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/sdc-pressure-outlet-checks/` contains checker, inputs, CSVs,
build/run/provenance/report and rejected-input/debug logs. Selected build and
CTest were not rerun. Prior SDC prescribed-pressure input-rejection evidence
is superseded. This establishes active routing and numerical inventories, not
large pressure-jump robustness, external solver parity or physical validation.
Inflow/walls and AMR/MPI SDC remain open.
