# Fixed-wall transport in coupled EB CTU — 2026-09-25

Implementation `2552b21`, based on `794b125` on `main`.
PeleC pin remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

The single-level EB transport flux/RHS/Euler/SSPRK2 paths now apply explicit
`slip_wall`, `no_slip_wall` and `isothermal_wall` policies. Walls have zero mass
and species exchange. Slip walls have zero tangential traction; stationary
adiabatic walls have zero energy flux. Tangential velocity derivatives respect
fixed wall velocity. Isothermal walls reuse the existing regular
`wall_heat_transfer_3d` kernel, with the heat sign oriented into the fluid and
wall-temperature thermodynamic validation. Apertures enter the divergence and
boundary ledger once. The full gas driver forwards wall temperature with built-in
stage donor preparation, retaining atomic state/ledger publication.

Custom transport providers remain restricted to inflow/outflow in the full driver.
Periodic transport seams, moving walls, wall timestep/input/restart application
wiring and full reference parity remain pending. Physical wall applications are
not enabled by this library change. The caller must supply extended EB metrics
consistent with wall-reflected donors, and an appropriate stable interval.

## Focused evidence

GNU Fortran 13.3.0, Python 3.12.3, existing Ninja Debug/OpenMPI configuration in
`build/spray-amr-reactive-debug`. No compiler/backend changes.

- Existing `unit_reactive_eb_transport_3d` passed in 0.97 s.
- Final `unit_eb_transport_halo_3d` passed in 15.86 s (15.91 s CTest total).
  It adds 33 open-face comparisons across three wall labels and all directions,
  using full-fluid and planar cut-cell geometries. Uniform moving-state face
  fluxes agree with the regular Fortran transport operator at `3e-12` scaled
  tolerance; wall mass/species fluxes and stationary adiabatic work are zero.
  Hot-wall flux points into the fluid. Six coupled reacting CTU cases close the
  boundary/wall/chemistry inventory at `2e-12` scaled tolerance. Invalid wall
  temperature retains the input state/T and caller boundary ledger.

The first combined test run had one pass and one failure: the cut-cell coupled
wall fixture used outward-extrapolated geometry while state ghosts reflected at
all six walls. On this four-cell domain, six halo layers reflect repeatedly and
can reach a covered donor from an open extrapolated metric cell. The existing
hydro guard correctly rejected it. The fixture now builds reflected nodal
geometry with an unchanged physical center; the final test passed. No numerical
bounds, tolerances or guards were loosened. Compiler warnings remain, including
the existing executable-stack warning for the test's internal provider.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_transport_halo_3d test_reactive_eb_transport_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_(eb_transport_halo_3d|reactive_eb_transport_3d)$' --output-on-failure
# After correcting only the wall fixture:
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_transport_halo_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_transport_halo_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, contract and whitespace checks passed. Mechanism SHA unchanged:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Artifacts: `build/eb-wall-transport-checks/` containing `build*.log`,
`LastTest-initial.log`, `diagnosis.log` and `LastTest-final.log`.
The temporary fixture diagnostic print was removed before the final build.
No application, full matrix, installed-binary audit, external-PeleC comparison or
physical wall-flow validation was performed here.
