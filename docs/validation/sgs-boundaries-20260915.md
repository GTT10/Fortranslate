# Physical-boundary SGS — 2026-09-15

Source: parent `23e0684` plus this commit's LES implementation, spray driver,
tests/checker and docs. GNU Fortran 13.3 Debug in `build/port-spray-debug`.

The existing constant-coefficient deviatoric Smagorinsky flux now supports
periodic, outflow, slip and stationary adiabatic no-slip faces. Two ghost layers
supply the transverse gradients. Only physical faces and internal faces enter
the flux loop and maximum-diffusivity estimate; artificial outer ghost faces
are excluded. A face exchange ledger is integrated through SSPRK2 and both LES
split half-steps. The periodic calculation keeps its original face loop.
Failed updates preserve gas state, temperature and the accepted surface ledger.

## Focused results

Five tests passed in 6.31 seconds:

- `spray_les_boundary`: uniform U=20 m/s flow between stationary y walls gives
  modeled wall viscosity `rho (Cs dx)^2 |2U/dy|` and net traction
  `-4 mu_t U A/dy`. The numerical boundary rate matches within relative 1e-12;
  stationary walls perform no work. The volume divergence matches the surface
  rate, RK momentum exchange agrees within relative 1e-10, and an excessive
  timestep rolls back state/temperature with zero committed ledger.
- `spray_sgs_boundary_restart`: Cs=0.17 with 10 m/s sinusoidal gas shear, inert
  water/drag and both wall/outflow cases. Stop after two steps and resume to
  0.6 microseconds; gas/parcel CSVs and final history rows match exactly.
  Mass/momentum/energy residual bounds and changed-face restart rejection pass.
- Existing `spray_les_grid`, `spray_les_zero` and `spray_les_reject`: periodic
  conservation/dissipation, exact zero-coefficient identity and rollback pass.

Artifacts are under `build/port-spray-debug/tests/spray_sgs_boundaries/` and
`build/port-spray-debug/Testing/Temporary/LastTest.log`.

```sh
cmake --build build/port-spray-debug --target test_spray_core pelef_spray_3d --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(les_(boundary|grid|zero|reject)|sgs_boundary_restart)$' --output-on-failure
python3 -m py_compile tools/check_spray_boundaries.py
python3 tools/check_project_contract.py --source .
python3 tools/check_mechanism_provenance.py --bundle build/port-methanol-reference/ffcm1_red.json --source build/port-methanol-reference/ffcm1_red.yaml
```

Project contract, pinned source provenance, Python syntax and diff whitespace
checks pass. No full matrix, experimental wall-flow benchmark or reacting
physical-boundary spray validation was run. Wall damping, wall functions,
resolved near-wall qualification and general geometry remain open.
