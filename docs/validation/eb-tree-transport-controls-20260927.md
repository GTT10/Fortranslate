# EB tree transport controls — 2026-09-27

Implementation: `46db17e1`, based on `980b20aa`.

The shared EB tree application previously hard-coded all four transport
operators on whenever gas transport was enabled. Its namelist now exposes
viscosity, thermal conduction, species diffusion and barodiffusion controls,
with the existing all-on defaults. They reach split and SDC, serial/MPI gas
and MPI spray, common/recursive stepping and automatic timestep estimates.
Barodiffusion is disabled when species diffusion is disabled. Isothermal-wall
input rejects disabled thermal conduction. Existing checkpoint configuration
fields bind these settings; the record format is unchanged.

## Completed focused checks

GNU Fortran 13.3.0 Debug, OpenMPI 4.1.6, MPI/spray ON,
tests/SUNDIALS OFF, selected pinned heptane Simple bundle:

```sh
python3 tools/check_mechanism_provenance.py --bundle build/heptane-3sp-simple-reference/heptane_3sp.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_reactive_eb_tree_3d_selected pelef_mpi_reactive_eb_tree_3d_selected --parallel 2
python3 build/run_spray_transport_isolation.py
python3 build/check_spray_transport_controls.py
python3 -m py_compile tools/check_project_contract.py build/run_spray_transport_isolation.py build/check_spray_transport_controls.py
python3 tools/check_project_contract.py --source .
git diff --check
```

All three targets build. Logs: `build/spray-transport-controls-build.log` and
`build/spray-transport-controls-gas-build.log`. Existing real-comparison,
array-temporary and executable-stack warnings remain.

Four two-step 4³ single-rank spray runs complete: viscosity only, conduction
only, species diffusion only, and the unmodified default input. Default gas
and particle CSVs are byte-identical to `build/spray-face-transport/`.
Artifacts and exact command arrays: `build/spray-transport-isolation2/`.

A conduction-only checkpoint after step one resumes on two ranks with gas
CSV byte-identical to the continuous single-rank run. Changing viscosity on
restart fails with `Gas snapshot restore failed`, exit 1, without producing
a final gas CSV. Artifacts: `build/spray-transport-controls-check/`, including
`commands.json` and `result.json`. This verifies effective transport settings
are bound by the existing gas checkpoint contract. Particle restart equality
was not independently compared in this increment.

The initial pre-implementation viscosity-only attempt rejected its unknown
namelist controls (`build/spray-transport-isolation/viscosity/fortran.log`).
It is not counted as a completed run.

## Reference diagnostic limit

Three C++ runs also complete with corresponding `diffuse_*` switches, but
their comparisons are **not operator-equivalent parity tests**. Pinned
`Source/Diffusion.cpp` computes all diffusion and then zeroes output flux
components. Either `diffuse_temp` or `diffuse_enth` keeps the complete energy
flux, including conduction, viscous work and species enthalpy transport.
Fortran's physical-operator controls instead include energy contributions
only from enabled operators. The exploratory comparison JSONs are retained
but do not establish a cause for the remaining full-transport residual.

No full configuration matrix, physical validation, fixed-mechanism rebuild,
recursive runtime or serial runtime was run in this increment. The selected
serial and MPI gas builds verify those compilation branches only.
