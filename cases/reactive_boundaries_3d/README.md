# Regular 3D physical boundaries

`pelef_reactive_3d` and the selected-mechanism front end accept six face
conditions in `boundary_faces`, ordered **x-low, x-high, y-low, y-high,
z-low, z-high**. Empty entries inherit `boundary_condition` (default
`periodic`). Opposing periodic faces must be paired.

- `inflow`: prescribed, time-independent reservoir state, specified by
  `inflow_temperature(6)`, `inflow_pressure(6)` and `inflow_velocity(3,6)`.
  Values use SI units and the same six-face order. The shared inflow composition
  defaults to the initial mixture; override it with `inflow_composition_count`,
  `inflow_composition_species` and `inflow_composition_mole_fractions` (species
  names, not mechanism indices). Both fixed and selected mechanisms resolve
  those names against their loaded species. For example:

  ```fortran
  boundary_faces='inflow','outflow','periodic','periodic','periodic','periodic'
  inflow_temperature(1)=1200, inflow_pressure(1)=101325
  inflow_velocity(:,1)=100,0,0
  inflow_composition_count=3
  inflow_composition_species='H2','O2','N2'
  inflow_composition_mole_fractions=.2,.2,.6
  ```

  This prescribes ghost states for the Riemann and transport fluxes; it is not
  a characteristic subsonic inlet or a prescribed numerical mass-flow rate.
  At intersecting inflow faces, z takes precedence over y over x for corner
  ghosts; wall reflections then apply. Avoid conflicting reservoir states at
  adjacent inlets unless this corner convention is intended. The timestep
  also accounts for inlet wave speeds and inlet transport diffusivity.
  To use the pinned nozzle-style Mach inlet instead, set `inflow_mach(face)`
  between 0 and 100 (default -1 retains the reservoir). The velocity points
  inward along the face normal with magnitude Mach times the inlet sound speed;
  `inflow_velocity(:,face)` must remain zero. Temperature is prescribed by
  `inflow_temperature`; composition comes from each stage's interior donor,
  including when a shared reservoir composition is supplied for other faces.
  `inflow_pressure(face)=0` extrapolates interior pressure; a positive value
  prescribes pressure (default 101325 Pa, maximum 1e15 Pa). For example:

  ```fortran
  inflow_mach(1)=0.2
  inflow_temperature(1)=1200
  inflow_pressure(1)=0
  ```

  Both regular fixed/selected front ends connect this mode through reaction,
  molecular transport and hydro. Before each step, padded current inlet states
  contribute to the hydro and diffusion timestep limits. This adds padded-grid
  work only when a Mach inlet is active. `mach_inlet.nml` combines extrapolated
  and prescribed-pressure Mach faces with chemistry and molecular transport.
  Other applications reject active Mach inputs until their dispatch is connected.
  This is a prescribed Mach ghost state, not GC-NSCBC or nozzle validation.
- `outflow`: copy the nearest physical cell into ghost cells (zero-gradient
  extrapolation). Optional `outlet_pressure(6)` entries select characteristic
  pressure extrapolation on these faces: zero keeps the default, positive
  values specify pressure in Pa (finite, at most 1e15). For example,
  `outlet_pressure(2)=102325` sets the x-high pressure target. Hydro and
  molecular transport share both ghost layers; composition is extrapolated
  from the interior. Invalid extrapolated states reject the step. This does
  not prescribe a backflow reservoir or establish nonreflection.
  `pressure_outlet.nml` exercises the full chemistry/transport/hydro path.
  Pressure input is enabled in regular, serial/MPI two-level AMR and serial/MPI
  patch-tree fixed/selected front ends.
- `slip_wall`: mirror the normal momentum and extrapolate thermodynamic and
  tangential states by reflection. The inviscid wall is impermeable and
  adiabatic; it exerts pressure forces on the gas.
- `no_slip_wall`: stationary, impermeable, adiabatic wall. The transport
  operator reflects all three velocity components, imposing zero face velocity
  while allowing viscous wall shear. The Euler operator reflects only normal
  momentum so its numerical dissipation does not add artificial wall friction.
- `isothermal_wall`: stationary no-slip wall with a prescribed
  `wall_temperature(6)` in kelvin (default 300). Molecular transport and heat
  conduction must be enabled. The temperature must lie in the loaded mechanism's
  thermodynamic range. Wall conductivity is evaluated at the prescribed wall
  temperature and adjacent-cell composition/pressure. The heat input is
  `k_wall*(T_wall-T_cell)/(dx_normal/2)`; only the physical face energy flux is
  changed. Mass/species fluxes remain zero and the wall shear is retained.
  This direct face condition does not extrapolate a possibly negative ghost
  temperature. The thermal timestep cap uses the adjacent cell's volumetric
  heat capacity and sums conductances at wall corners; half the normal
  interior diffusion allowance is reserved when these walls are present.
- `periodic`: wrap cell states from the opposite side.

PCM and characteristic PLM use the existing Rusanov, HLLC or PeleC Riemann
solver. Two ghost layers are filled at each SSPRK2 stage. Only the physical
cells are advanced. Stage-weighted physical surface fluxes are accumulated
in `Integrated boundary change`, ordered as density, three momenta, total
energy, then species densities, integrated over space and time. The
application subtracts these exchanges before checking conservation. Closed wall domains conserve mass. Stationary adiabatic walls conserve total
gas energy; moving or isothermal walls exchange energy through work or heat. Chemistry and molecular transport are included in the existing transactional
split step. A shared ghost-fill routine supplies viscous stresses, heat conduction,
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

`no_slip.nml` advances initially uniform tangential flow between stationary
y walls with periodic x/z faces. Viscosity transfers momentum to the walls
and converts kinetic energy to internal energy; total gas energy is conserved.

`inflow.nml` exercises reservoir injection with molecular transport and an
independent inlet composition. The boundary ledger accounts for both convective
and diffusive exchange with the reservoir.

`isothermal.nml` uses 350 K and 1100 K stationary y walls around an initially
1000 K gas. The boundary ledger includes heating/cooling, so total gas energy
can change consistently with the specified wall heat flux.

`wall_velocity(3,6)` prescribes constant tangential velocity in m/s for
`no_slip_wall` and `isothermal_wall` faces (zero by default). Molecular viscosity
must be enabled. For example, `wall_velocity(:,4)=250,0,0` moves the top y wall
in the x direction. The transport ghost velocity is reflected about the wall
velocity, with total energy adjusted to preserve internal energy. The energy
flux includes wall traction dotted with wall velocity, plus prescribed heat
conduction for an isothermal wall. Hydro remains impermeable. Normal wall motion
and nonzero velocity on other face types are rejected. At wall intersections,
ghost reflections are composed in x/y/z order.

`moving_wall.nml` exercises this work and momentum transfer between periodic
x/z faces and adiabatic y walls. See [focused evidence](../../docs/validation/moving-walls-20260915.md).

Moving geometry, characteristic inflow, general AMR/EB boundaries and particle boundary
coupling remain open. The padded-grid implementation prioritizes reuse of the
existing numerical fluxes; memory/scaling qualification remains open.

## External forcing

Regular fixed/selected gas applications accept `external_forcing(3)` in N/m³
(force per volume), default zero. It requires `characteristic_ppm` and `pelec`.
The momentum source equals this vector; energy uses the average of old/new
velocity dot force. `Integrated source change` reports volume sources separately
from `Integrated boundary change`; conservation diagnostics subtract both.

Run `pelef_reactive_3d cases/reactive_boundaries_3d/external_forcing.nml` for a
small uniform periodic example. The selected executable needs its usual named
composition and `thermo_model='selected'`. Nonzero forcing is rejected by
application families whose input reader has not enabled this feature.
These regular gas executables have no checkpoint interface. This option does
not establish full SDC integration or physical forcing validation.

## SDC

`time_integrator='sdc'` selects iterative source-coupled CTU, with
`sdc_iterations=2` by default. The default integrator remains `split`.
The current application path requires PPM/PeleC and periodic, inflow or
outflow boundaries. Molecular transport is supported at physical boundaries;
`outlet_pressure(6)` also supports prescribed pressure on outflow faces.
Inflow accepts the existing specified-state and Mach-number input fields.
Walls are not yet connected to SDC. Boundary flux
inventories are recorded separately from volume sources.
With `chemistry_enabled=.true.`, it uses forced adaptive implicit chemistry
and retains reaction sources between steps. The chemistry tolerances apply;
the absolute tolerance is scaled by the cell's starting density for the
species-density reactor. Molecular transport is evaluated as an old/new
SDC source when `transport_enabled=.true.`, using the existing transport flags
and stable timestep. It supports the external force input above. Unsupported combinations fail before stepping.

Run `pelef_reactive_3d cases/reactive_boundaries_3d/external_forcing_sdc.nml`
for the two-iteration example. Fixed/selected gas variants share the driver;
SDC is not yet enabled for AMR, MPI or EB. `sdc_reacting.nml` provides a
uniform reacting example with transport and external forcing; this is a
bounded numerical example, not physical ignition validation.
`sdc_outflow.nml` exercises a reacting hotspot with external force and outflow.
Set `transport_enabled=.true.` in that example to include molecular transport.
