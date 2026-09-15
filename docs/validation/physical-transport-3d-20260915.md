# Molecular transport at regular 3D physical boundaries

Parent: `f3f5632`. Physical outflow/slip-wall boundaries now support viscosity,
heat conduction, species diffusion and barodiffusion in the fixed/selected
regular 3D applications. See the [case guide](../../cases/reactive_boundaries_3d/README.md).

Two shared ghost layers supply the existing transport flux and species donor
limiter. Physical cells alone enter the Euler update; overflow checks and
EOS rejection are retained. Both SSPRK2 stages contribute surface budgets,
and the split driver accumulates both transport half-steps plus hydro fluxes.
A rejected step commits neither gas state nor its boundary ledger. Periodic
calls retain their original flux/update dispatch.

## Focused local evidence

GNU Fortran 13.3.0 Debug, `build/port-spray-debug`:

```sh
cmake --build build/port-spray-debug --parallel 8 --target \
  test_reactive_transport_3d test_reactive_boundaries_3d \
  test_reactive_boundary_transport_3d pelef_reactive_3d
ctest --test-dir build/port-spray-debug --output-on-failure \
  -R '^unit_reactive_(transport_3d|boundary_transport_3d|boundaries_3d)$'
ctest --test-dir build/port-spray-debug --output-on-failure \
  -R '^reactive_boundary_3d_transport$'
```

- Existing periodic transport and physical hydro tests passed.
- New transport test passed (1.68 s final run): each physical process separately
  and combined, wall impermeability/adiabaticity, open viscous work, surface
  conservation, species closure/positivity, an active donor limiter, rollback.
- The mixed-boundary application passed (0.81 s), exercising transport/hydro
  splitting and the accumulated conservation ledger.
- An initial new-test assertion did not detect enough net open viscous work
  with closed transverse faces. The final fixture opens x and y and keeps z
  periodic, producing nontrivial boundary work at the unchanged threshold.
  No numerical acceptance tolerance was loosened.

The existing FFCM1_Red Release configuration also built
`pelef_reactive_3d_selected` and completed `physical-transport-run/transport.nml`:
CH3OH/O2/N2 mole fractions 0.05/0.20/0.75, otherwise the committed transport
case. Its boundary-accounted normalized Euler residual was `6.305e-17`.
Chemistry was disabled in this selected smoke; the test is not fuel validation.

Provenance and project-contract checks and `git diff --check` passed. No Python
source was changed. Local build/test logs are under `build/port-spray-debug/`
with prefix `physical-transport-`. The selected build log and run directory
are under `build/port-spray-cvode/`. This is focused local evidence, not a full
configuration matrix, external field parity or project-completion claim.

Prescribed inflow, no-slip/isothermal walls, general AMR/EB physical boundaries,
particle coupling and the other roadmap work remain open.
