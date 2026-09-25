# Regular PCM/PLM explicit artificial viscosity — 2026-09-25

## Source and scope

Code `b3cf4bd`, parent `90cbad2`. Regular serial boundary hydro applies the
existing artificial-viscosity kernel to PCM/PLM fluxes at both SSPRK2 stages.
The full-step dispatcher routes periodic PCM/PLM with a positive explicit
coefficient through this path. Published mean face fluxes and boundary ledgers
include the correction. No-slip/isothermal faces suppress tangential momentum
corrections through the existing kernel. Failure retains the caller state.

Reference inspection used pinned PeleC
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8:Source/Hydro.cpp`:
regular `pc_adjust_fluxes` applies artificial viscosity before species-flux
normalization. This increment connects the existing Fortran kernel; it does
not establish complete PeleC time-integration parity.

Omitted and zero PCM/PLM coefficients retain the previous zero policy. Explicit
non-PPM namelist input remains unavailable. Spray wrappers, MPI/AMR propagation,
consistent default/input/restart policy and EB geometry remain future work.
The full-step test disables reactions and molecular transport.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug, Python 3.12.3; existing tests-enabled build
`build/spray-amr-reactive-debug`. Focused CTest **1/1 passed**, 2.37 s
(2.41 s elapsed). No failures or skips. Compiler real-comparison warnings remain.

The 4-cubed fixed seven-species fixture exercises PCM and characteristic PLM,
periodic boundaries and a no-slip/outflow pair. It checks omitted/zero exact
identity, nonzero response to coefficient 0.037, reconstruction of the update
from published mean fluxes, boundary inventory closure, species closure,
full-step dispatch exact parity, and negative-coefficient rollback. Numerical
closure bounds are 2e-13 with the scales recorded in the test. This is focused
library evidence, not external solver parity or physical validation.

Mechanism provenance, Python syntax, project contract and `git diff --check`
passed. Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_pcm_plm_viscosity_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_pcm_plm_viscosity_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/pcm-plm-viscosity-checks/build.log` and copied `LastTest.log`.
