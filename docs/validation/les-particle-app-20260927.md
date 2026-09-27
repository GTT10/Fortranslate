# Particle LES SDC application — 2026-09-27

Implementation `ea6c7ede`, based on `96f51cda`, local main, version 0.246.0.

## Change

The MPI EB spray frontend reads optional `eb_spray_les` controls from its spray
input. Enabling LES requires SDC. Smagorinsky, WALE, Vreman and Dynamic settings
feed the regular-geometry LES adapter through the common/recursive particle
wrappers, including their existing wall-film branches. The adapter rejects
non-regular geometry. A pre-step LES evaluation validates the current state
and limits the timestep using `les_cfl`, the largest diffusivity and the current
smallest patch spacing. The bound is conservative for recursive stepping.

Enabled LES settings are serialized into the existing bound spray/gas/particle
checkpoint context. Changed or removed LES settings reject continuation.
Disabled LES adds no context record and uses the existing callback-free path.
The application documentation gives the namelist and geometry restrictions.

## Independently completed checks

Artifacts: `build/les-particle-app-checks/` (`provenance.log`, `build.log`,
`build-final.log`, `ctest.log`, `LastTest.log`). Per-run input/log/checkpoint/CSV
artifacts:

- `build/spray-amr-reactive-debug/les-particle-app-common/les_q_8yrw70/`
- `build/spray-amr-reactive-debug/les-particle-app-recursive/les_yzkrt30u/`

Existing fixed-mechanism Ninja Debug build, GNU Fortran 13.3.0, Open MPI 4.1.6
(`/usr/bin/mpiexec`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_les_sdc_app_(common|recursive)$' --output-on-failure
python3 -m py_compile tests/application/check_les_particle_sdc_app.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, no failures/skips: common 24.34 s, recursive 52.65 s,
CTest total 77.05 s. Build, mechanism provenance, Python syntax, project
contract and whitespace checks passed. The initial compile rejected long
argument lines; wrapping them resolved the failure. The executable-stack
linker warning remains. Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

The fixed-mechanism two-level fixture has fully regular nodal geometry, outflow
faces, molecular transport, two water parcels and disabled chemistry. Enabled
Smagorinsky LES changes the gas CSV relative to disabled LES. Combined inventory
residuals are 1.7604525331008246e-15 (common) and 1.2477363672874494e-15
(recursive), below the 3e-12 normalized criterion. A checkpoint saved after
step 1 on 2 ranks continues on 4 ranks with identical final gas/budget bytes
and particle rows to the continuous 4-rank case. Changed `les_cs`, removed LES,
and valid cut-cell input geometry are rejected with the expected diagnostics.

## Limits

The selected executable, other LES model inputs, a timestep actually shortened
by the diffusion bound, active chemistry, films, regridding, other boundary
types, installed binaries, full configuration matrix, external PeleC parity
and physical validation were not tested here. All four closures have separate
particle-library evidence. This remains a regular-cell LES model, not a
cut-cell LES extension.
