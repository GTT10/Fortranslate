# Spray gas reservoir inlet — 2026-09-15

Source: parent `9449f1f` plus this commit's spray/LES inlet plumbing, input
fixture, checker and documentation. GNU Fortran 13.3 Debug in
`build/port-spray-debug`, with spray enabled.

## Implemented scope

Six-face reservoir inflow supplies named mixture composition, temperature,
pressure and velocity to hydro, molecular transport and SGS ghost states.
The inlet mixture defaults to the initial gas mixture; named overrides resolve
against the loaded fixed or selected mechanism. Active inlet wave speed and
molecular diffusivity cap the application timestep. Source parameters and the
resolved mixture are included in checkpoint context; existing physical-boundary
schemas already contain the complete surface inventory.

A parcel that leaves through an inlet is treated as particle outflow. Gas
inflow does not manufacture parcels or alter the internal injection schedule.
This is a prescribed reservoir boundary, not a characteristic inlet or a
numerically fixed mass-flow boundary.

## Focused evidence

`spray_inflow_restart`, existing `spray_sgs_boundary_restart` and
`spray_application_cli` passed 3/3 in 23.34 seconds.

The inlet fixture supplies 1100 K, 101325 Pa, 200 m/s H2/O2/N2 gas into an
initially pure-N2 domain containing two evaporating water parcels. Molecular
transport and Cs=0.17 are enabled; chemistry is disabled. Stop after two steps
and restart to 0.6 microseconds:

- Gas and parcel CSVs and the final history row match uninterrupted integration
  exactly.
- Net boundary mass is `3.24085191656049596e-10` kg. Mass and energy residuals
  are `2.15583947146632838e-23` kg and `-4.87890977618476995e-18` J.
- H2 enters the gas field, establishing that the independent named inlet
  composition reaches the solution. All three momentum residuals satisfy the
  checker's fixed bounds.
- Changing inlet temperature to 1101 K rejects restart before creating outputs.

Artifacts:
`build/port-spray-debug/tests/spray_inflow_restart/spray-boundaries-1i8xnkik/`.
Project contract, pinned FFCM1_Red source provenance, Python syntax and
`git diff --check` pass.

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(inflow_restart|application_cli|sgs_boundary_restart)$' --output-on-failure
python3 -m py_compile tools/check_spray_boundaries.py
python3 tools/check_project_contract.py --source .
python3 tools/check_mechanism_provenance.py --bundle build/port-methanol-reference/ffcm1_red.json --source build/port-methanol-reference/ffcm1_red.yaml
```

No full compiler/backend matrix, reacting selected-mechanism inlet benchmark,
resolved nozzle comparison or experimental spray validation was run. Boundary
particle generation, characteristic inflow and moving/isothermal spray walls
remain open.
