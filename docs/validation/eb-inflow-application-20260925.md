# EB CTU inflow application input and restart — 2026-09-25

Implementation `82e5675`, based on `bf89878` on `main`.
PeleC pin: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

Fixed and configure-time selected single-level gas front ends now accept six
explicit inflow/outflow labels and per-face density, pressure and Cartesian
velocity. Inlets use the initial mixture composition. Legacy labels remain the
default and retain old hydro/transport behavior. Mixed labels and unsupported
methods/policies reject; front ends without explicit reader opt-in reject the new
boundary mode, preventing silent use by spray or tree callers.

The shared application constructs thermodynamically checked prescribed states,
forwards them through coupled CTU, and rebuilds one-layer donors for each timestep
estimate. Only open-face donors contribute signal-speed and diffusivity bounds.
Checkpoint solver identity appends `boundary-v1`, six labels and all inlet values;
legacy identity text remains unchanged. The identity buffer is enlarged to 4096
characters to hold these fields. An input example is
`cases/reactive_eb_3d/plm_ctu_inflow.nml`.

## Evidence

GNU Fortran 13.3.0, Python 3.12.3, Ninja Debug, existing OpenMPI wrappers.
Fixed tests-on build: `build/spray-amr-reactive-debug`.
Selected tests-off build: `build/spray-amr-reactive-selected-debug`,
`mechanisms/h2o2_full.json`. No install prefixes or backend settings changed.

- Config/checkpoint unit tests: **2/2 passed**, 0.24 s total (0.06 s config,
  0.14 s checkpoint), no failures or skips.
- Fixed application checker: **5 successful processes, 6 expected rejections**.
- Selected H2O2 application checker: **5 successful processes, 6 expected rejections**.
- Each checker covers inert, reacting, coupled, checkpoint-stop and resume runs;
  continuous and resumed final CSVs match byte-for-byte, component inventories
  satisfy the existing checks, and boundary-label or inlet-velocity changes reject
  at restart without creating the requested CSV. Method/solver/invalid-inlet
  combinations reject during input validation.
- One additional fixed first-step run changes y-inlet velocity from 200 to
  400 m/s. Its timestep decreases from `1.0028096763855324e-10` to
  `9.2833426932924000e-11` s. It uses a separate `faster.chk` and output file.
  This directly checks the inlet acoustic restriction; it does not separately
  isolate a diffusion-dominated timestep.

Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_eb_3d test_reactive_eb_3d_config test_reactive_eb_3d_checkpoint --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_eb_3d_(config|checkpoint)$' --output-on-failure
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 python3 tools/check_eb_ctu_application.py --executable build/spray-amr-reactive-debug/pelef_reactive_eb_3d --work-dir build/eb-inflow-application-checks/fixed --inflow
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 python3 tools/check_eb_ctu_application.py --executable build/spray-amr-reactive-selected-debug/pelef_reactive_eb_3d_selected --work-dir build/eb-inflow-application-checks/selected --inflow --selected
python3 -m py_compile tools/check_eb_ctu_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, Python syntax, project contract and whitespace checks passed.
Compiler warnings remain; no build/runtime failures occurred.
Mechanism SHA remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Artifacts under `build/eb-inflow-application-checks/`:
`{build.log,build-selected.log,LastTest.log}`, fixed `run-lqpj0d4r`, selected
`run-z3k25yl0`. The additional first-step input/output/log are `faster.*` in the
fixed run directory, derived from `stopped.nml` with the velocity and artifact
names changed. No broader matrix, installed-binary audit, external-solver parity,
or physical inlet-flow validation was performed.

Periodic/wall transport, independent inlet mixture selection, spray/tree/MPI
boundary input and complete reference parity remain separate work.
