# Regular 3D physical boundaries

`pelef_reactive_3d` and the selected-mechanism front end accept six face
conditions in `boundary_faces`, ordered **x-low, x-high, y-low, y-high,
z-low, z-high**. Empty entries inherit `boundary_condition` (default
`periodic`). Opposing periodic faces must be paired.

- `outflow`: copy the nearest physical cell into ghost cells (zero-gradient
  extrapolation). This is not characteristic or non-reflecting outflow and
  does not prescribe a backflow state.
- `slip_wall`: mirror the normal momentum and extrapolate thermodynamic and
  tangential states by reflection. The inviscid wall is impermeable and
  adiabatic; it exerts pressure forces on the gas.
- `periodic`: wrap cell states from the opposite side.

PCM and characteristic PLM use the existing Rusanov, HLLC or PeleC Riemann
solver. Two ghost layers are filled at each SSPRK2 stage. Only the physical
cells are advanced. Stage-weighted physical surface fluxes are accumulated
in `Integrated boundary change`, ordered as density, three momenta, total
energy, then species densities, integrated over space and time. The
application subtracts these exchanges before checking conservation. Closed
wall domains conserve mass and total energy, but the wall can change gas
momentum. Chemistry and molecular transport are included in the existing transactional
split step. The same boundary ghosts feed viscous stresses, heat conduction,
species diffusion and barodiffusion. Surface budgets include the two RK stages
of both transport half-steps. Species donor limiting and enthalpy correction
use the existing conservative transport operator.

```sh
mkdir -p build/boundary-example
cd build/boundary-example
../port-spray-debug/pelef_reactive_3d ../../cases/reactive_boundaries_3d/outflow.nml
../port-spray-debug/pelef_reactive_3d ../../cases/reactive_boundaries_3d/closed_reacting.nml
```

`outflow.nml` combines open x faces, periodic y faces and closed z faces.
`closed_reacting.nml` advances a reacting hotspot in a closed box. They exercise
regular-grid hydro/chemistry and boundary-accounted conservation, not physical
validation against an experiment. The periodic entropy-wave exact-solution
error is not reported for physical-boundary runs.

`transport.nml` exercises all molecular transport processes with mixed faces.
The stationary slip wall has zero normal heat/species flux and zero tangential
traction; its normal viscous traction can change momentum. Outflow extrapolates
zero normal state gradients, while tangential derivatives still contribute to
stress and its energy flux.

Prescribed inflow, no-slip/isothermal walls,
characteristic boundaries, general AMR/EB boundaries and particle boundary
coupling remain open. The padded-grid implementation prioritizes reuse of the
existing numerical fluxes; memory/scaling qualification remains open.
