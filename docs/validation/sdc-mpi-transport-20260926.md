# MPI SDC molecular transport stage — 2026-09-26

Implementation `a62a754c`, base `55a41f3b`; PeleF 0.246.0.
Pinned PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

## Behavior

`compute_sdc_mpi_amr_tree_transport_3d` returns owner-local transport rates,
dt-integrated face fluxes and the collectively reduced minimum limiter theta.
Its four flags control viscosity, conduction, species diffusion and
barodiffusion. It uses existing MPI viscous ghost exchange with MC coarse
interpolation and existing pressure/inlet/thermal/moving-wall handling.

The read-only wrapper reuses the ordinary transport routine's control and
boundary validation plus face loop through a collectively agreed stage mode.
This mode skips initial restriction/EOS recovery, state advancement, second
Euler evaluation and final reflux. The wrapper holds private owner-local
copies; no complete global state is gathered. Outputs publish only after
collective success. Ordinary transport behavior is preserved.

## Completed focused checks

Artifacts: `build/sdc-mpi-transport-checks/` (local).
GNU Fortran 13.3.0, Ninja Debug, existing tests-enabled
`build/spray-amr-reactive-debug`; configured Open MPI 4.1.6 with absolute
`/usr/bin/mpiexec`. No compiler/backend changes.

- Provenance passed: `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`;
  SHA `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
- Build command: `env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build
  build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_transport_3d
  --parallel 2`. Initial build failed at the 132-column Fortran limit in the
  expanded flag array (`build.log`). Wrapping the expression resolved it;
  final build passed (`build-final.log`). Existing warnings remain.
- Focused CTests: **2/2 passed**, 11.43 s total (`ctest.log`).
  New `unit_mpi_amr_patch_tree_sdc_transport_3d_4`: 1.60 s.
  Existing `unit_mpi_amr_patch_tree_transport_3d_2`: 9.78 s.
- New test: four ranks with ownership `[0,1,0]` and two idle ranks;
  nonuniform three-level state, all four transport flags. Periodic and
  outflow/isothermal-moving-wall cases yield exact serial/MPI rates, all three
  integrated flux arrays and theta. State/temperature remain unchanged and
  unowned output fields remain absent. Invalid transport data on one owner
  rejects collectively and preserves rates, fluxes and theta.
- Command: `env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH
  HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest
  --test-dir build/spray-amr-reactive-debug
  -R '^unit_mpi_amr_patch_tree_(sdc_transport_3d_4|transport_3d_2)$' --output-on-failure`.
- Python syntax for `tools/check_project_contract.py`, project contract,
  Fortran line lengths and final `git diff --check` passed.

No coupled MPI SDC iteration/application, level-time interpolation, selected
build, restart/rank change, full matrix, installation audit, full-solver pinned
comparison or physical validation was performed. The existing ordinary test
checks the shared update path; it is not a coupled SDC solve.
