# EB face halos and accepted flux integrals — 2026-09-17

## Source and behavior

Implementation `571b193aec677362adfacbe2d0970a049a862216`, based on
`aa05a48` on local `main`.

`reactive_eb_riemann_fluxes_3d` aliases the existing outflow flux API, now with
optional `ghost_state` and `ghost_temperature`. Both must be supplied together,
with one halo cell on every spatial side (0:n+1). Only donors adjacent to open
boundary faces are read; closed-face, interior, edge and corner halo slots are
ignored. Donors enter the left/right Riemann states according to face orientation.
Rusanov, HLLC and PeleC solvers retain their existing thermodynamic checks.
Absent halos retain extrapolated-outflow behavior.

`advance_reactive_eb_euler_3d` accepts the same optional halos and optional
`integrated_flux` of type `eb_integrated_flux_3d`. The flux returned after a
successful update is dt times the actual Riemann flux, per open face area, with
the existing normal-index lower bounds. Any failed update leaves that flux
object unchanged and returns the original conserved state and temperature under
the existing valid-shape rollback convention.

The conservative RHS and EB pressure-wall source are unchanged. This supplies
the donor and accepted-flux interfaces needed by patch stepping; it does not
construct coarse-grid donors or integrate an AMR hierarchy. Halos and returned
fluxes are currently exposed by the plain Euler step. State/flux redistribution
wrappers retain their previous interfaces. The explicit small-cell timestep
restriction still applies to plain Euler. Transport halos, subcycling and MPI
EB stepping remain open.

## Independently completed checks

GNU Fortran 13.3.0 Debug, Open MPI 4.1.6, Ninja; matching existing build
`build/selected-tree-mpi-debug`, MPI enabled and spray disabled. Fixed core and
selected EB runtime compiled. Execution used fixed elementary chemistry; no MPI
or application matrix was run. No Python source changed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 2 --target test_eb_reactive_ghost_3d test_reactive_eb_transport_3d pelef_selected_reactive_eb_3d_runtime
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 2 --target test_eb_reactive_ghost_3d
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^(unit_eb_reactive_ghost_3d|unit_reactive_eb_transport_3d)$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Final **2/2 PASS, 0 failures, 0 skips, 1.11 s**:

- New halo unit: 0.20 s. Regular and cut-cell 4-cubed geometries with all three
  Riemann solvers. An outflow-equivalent halo yields bit-exact legacy fluxes;
  modified donors match independent directional Riemann calls on both sides of
  each axis. Unused and closed-face NaN halo values are ignored. Accepted Euler
  flux integrals match the flux used for the update, and global conserved changes
  equal outer-boundary flux plus EB wall contribution for every component.
  Invalid donors preserve state, temperature and previously returned flux;
  missing partner arrays and wrong halo dimensions fail with zero flux output.
- Existing coupled transport unit: 0.88 s. Established hydro/transport/reaction,
  geometry and rollback paths pass with absent optional halos.

Provenance and project contract passed. Mechanism SHA remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Logs: `build/eb-ghost-checks/{build,final-build,ctest}.log`.
An earlier pass preceded the added closed-face-NaN and shape checks; only the
final run is counted. This is focused evidence, not complete AMR integration,
convergence, external PeleC parity or physical validation.
