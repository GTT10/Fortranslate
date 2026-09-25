# Boundary CTU step with conserved sources — 2026-09-26

Implementation `dfcf6743`, parent `4eb6d956`, local main, version 0.246.0.

`advance_reactive_boundary_hydro_3d` accepts optional conserved source rates
with four caller-supplied halo layers, matching the existing CTU flux contract.
The option requires `characteristic_ppm` with `pelec`. The predictor receives
the complete source array; the final candidate adds `dt*source` only in physical
cells. Rates are held fixed during this step. No source boundary values are
inferred from gas boundary labels.

`boundary_change` continues to report boundary flux only. The inventory equation
is final minus initial = boundary change + integrated physical-cell source.
Validation or EOS failure preserves state and temperature and returns zero
boundary change and optional flux outputs.

## Verification

`unit_reactive_ppm_step_3d` passed **1/1**, 3.50 s test / 3.55 s total.
The new nonuniform outflow fixture checks boundary/source inventory separation
for every conserved component. An excessive negative energy source checks full
state/temperature rollback and zero boundary change. Existing periodic source,
hydro conservation and splitting checks also ran. No failures or skips.
Provenance, Python syntax, project contract and whitespace passed.
GNU Fortran 13.3.0, Python 3.12.3, existing tests-enabled Ninja Debug build;
existing compiler warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_step_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_step_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/ppm-boundary-source-step-checks/{provenance.log,build.log,ctest.log,unit.log}`.
The underlying source/flux oracle is recorded in
[the 38-case comparison](ppm-conserved-source-20260926.md).
This is a focused library update check. Application source generation, MPI
source exchange, time-dependent quadrature/SDC and physical validation remain
open.
