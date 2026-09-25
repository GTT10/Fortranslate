# AMR tree gas/spray coefficient input and restart — 2026-09-25

## Implementation

Code `c71c95f`, parent `63be6a1`. Fixed/selected serial and MPI tree gas front
ends and MPI AMR spray front ends enable `difmag` in the gas namelist. Shared
application drivers forward the resolved value to the coupled operators.
`resolved_reactive_difmag_3d` centralizes automatic selection for these drivers
and the shared tree checkpoint context builder: `-1` means 0.1 for PPM/pelec,
zero otherwise; explicit nonnegative values retain their value.

Existing input validation requires PPM for an explicit coefficient and rejects
nonfinite/invalid negative values. The checkpoint builder canonicalizes the PPM
coefficient before serializing the configuration. Equivalent automatic/explicit
settings share context, and changed values reject. Old PPM checkpoints storing
an unresolved -1 sentinel are incompatible. Default PCM/PLM context remains
unchanged. No gas or particle payload format change was introduced.

The separate older two-level AMR family still lacks explicit coefficient
propagation, input enabling and restart binding.

## Completed verification

GNU Fortran 13.3.0, OpenMPI 4.1.6, Ninja Debug, Python 3.12.3. Existing
fixed tests-enabled `build/spray-amr-reactive-debug`; selected gas execution uses
the generated elementary-mechanism fixture targets in that build, not a
separately qualified selected installation.

**2/2 CTests passed**, gas restart 43.45 s, spray restart 6.33 s, total elapsed
49.82 s. Build, provenance, Python syntax, project contract and whitespace
checks passed. No failed checks or skips in this increment. Compiler warnings
remain; the build's MPI runtime dependency audit also passed (ABI 40).

- Gas: fixed and selected mechanisms, recursive timesteps, 2-to-3-rank and
  serial/MPI interchange, one-rank continuation, exact final CSV/budget output,
  and byte-identical serial/MPI checkpoints. Automatic -1 resumes the explicit
  0.1 checkpoint with identical outputs in serial and MPI. Explicit 0.2 fails
  with the expected context diagnostic, preserving existing output. Existing
  subcycling, physics, truncation and ledger guards remain active.
- Spray: fixed mechanism, recursive timesteps, timed injection and regridding,
  2-to-4-rank continuation at the existing relative 2e-12 / absolute 1e-25
  snapshot tolerances. Automatic -1 and explicit 0.1 have exactly equal
  four-rank gas/particle records. Changing to 0.2 fails with the checkpoint
  diagnostic before producing final CSV. Existing context/truncation/overwrite
  guards also pass.

Selected spray input and the shared driver were changed, but a selected spray
executable was not built or run here. No full matrix, installed audit, external
PeleC parity, detailed-fuel or physical validation was performed.

## Commands and artifacts

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_amr_tree_reactive_3d pelef_mpi_amr_tree_reactive_3d test_selected_elementary_amr_tree_3d test_selected_elementary_mpi_amr_tree_3d pelef_mpi_amr_spray_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(regression_mpi_amr_tree_difmag_restart_3d|spray_mpi_amr_difmag_restart)$' --output-on-failure
python3 -m py_compile tools/check_mpi_amr_tree_restart_3d.py tests/application/check_mpi_amr_spray_restart.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Logs: `build/difmag-amr-input-checks/{build.log,test.log,LastTest.log}`.
Gas fixtures and per-invocation logs:
`build/spray-amr-reactive-debug/tests/mpi_amr_tree_difmag_restart_3d/{fixed,selected}/`.
Spray fixture:
`build/spray-amr-reactive-debug/tests/spray_mpi_amr_difmag_restart/restart_uzgdka3h/`.
