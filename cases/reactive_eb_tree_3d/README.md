# Multilevel EB reacting tree

Build and run the fixed H2/O2 application:

```sh
cmake --build build/debug --target pelef_reactive_eb_tree_3d
./build/debug/pelef_reactive_eb_tree_3d cases/reactive_eb_tree_3d/input.nml
```

The supplied case has a 4-cubed root, two children, and grandchildren with
separate parents (three levels, ratios two). A planar wall cuts z=0.31.
Two 1e-9 steps run reaction, molecular transport and HLLC hydro.

## Artificial viscosity

Gas executables accept finite nonnegative `difmag` in `&eb_tree` (for example,
`difmag=0.037`). Default 0.1 matches the reference coefficient; explicit zero retains the previous policy; -1 is invalid.
The coefficient reaches common and recursive hydro and must match on restart.
MPI EB spray also forwards the coefficient through its gas update; see
[MPI spray evidence](../../docs/validation/mpi-eb-spray-difmag-20260925.md). See
[input/restart evidence](../../docs/validation/eb-tree-difmag-input-20260925.md)
for the checked configuration and remaining limitations.

## MPI execution

### Mach inlet

Fixed/selected serial and MPI gas executables and MPI EB spray accept six-face inlet arrays
in `&eb_tree`, ordered low-x, high-x, low-y, high-y, low-z, high-z. For example:

```fortran
 inflow_mach(1)=0.2,
 inflow_temperature(1)=400,
 inflow_pressure(1)=120000,
```

Mach 0..100 enables an inward normal inlet with the adjacent gas composition.
Temperature must be positive and EOS-valid. Pressure zero (default) extrapolates
the adjacent pressure; a positive value prescribes it. Negative Mach (default
-1) keeps outflow; temperature and pressure must remain zero on inactive faces.
All values must be finite. Common and recursive gas stepping, molecular boundary
fluxes and automatic timestep selection use the controls. Active inlet controls
are saved and compared during restart, including rank-changing continuation.
MPI EB spray forwards the same controls to its gas update; particle boundary
behavior remains governed by its separate particle settings. See
[gas application evidence](../../docs/validation/eb-mach-application-20260925.md)
and [spray coupling evidence](../../docs/validation/eb-mach-spray-20260925.md).

### Launch

With `PELEF_ENABLE_MPI=ON`, the fixed-mechanism executable
`pelef_mpi_reactive_eb_tree_3d` accepts the same `input.nml` and output/checkpoint
paths:

```sh
cmake --build build/mpi-debug --target pelef_mpi_reactive_eb_tree_3d
mpiexec -n 2 ./build/mpi-debug/pelef_mpi_reactive_eb_tree_3d cases/reactive_eb_tree_3d/input.nml
```

This entry point accepts static, prescribed scheduled or temperature-adapted
topology, fixed or automatic root intervals and explicit/implicit chemistry. Automatic estimates use collective
remote-donor maxima and the minimum patch hydro/diffusion bound, including hidden
parent cells. Scheduled regridding uses sparse requests to old owners, conservative
coarsening/PCM refinement and new whole-patch ownership. It preserves the root
lattice and immutable geometry, as in serial. `scheduled.nml` exercises full
coarsening and renewed refinement. `adaptive.nml` uses distributed exposed-fluid
temperature votes and deterministic tiling with solid-only patches removed.
Completed adaptive layouts are restored from the same V4 checkpoint envelope.
`subcycling=.true.` selects recursive hydro/transport and scales automatic
root intervals by cumulative refinement ratios (squared for diffusion).
Checkpoint stepping mode must match on restart. Each rank holds only its assigned patch
U/T during initialization and stepping; geometry is replicated. Initial children
fetch only their parent sections and use the serial PCM kernel on their owners.
The root patch owner constructs its seed U/T during each initialization or
restart topology setup and releases it when initialization returns. Other ranks
pass empty seed arrays. MPI checkpoints use a V2 byte-length index; on restart,
owners seek directly to their patch records and restore only local U/T.
Legacy V1 files from serial saves remain readable through a bounded-memory root scan.
State CSV writes use owner-local MPI-IO chunks in serial row order;
the small budget CSV is written by rank zero. Checkpoint writes stage existing
serial patch records on their owners and publish them directly with MPI-IO.
Checkpoints can continue with a different rank count or the serial
executable. State CSVs are deterministic in the focused fixture; budget ledgers
can differ in roundoff when MPI reduction grouping changes. Serial restart accepts
both formats and checks indexed record lengths. Scaling qualification remains future work.

The V2 tree record keeps the V1 topology and inventory header, followed by one
positive byte length per node, in node order, before the first patch record.
Writers use a 20-column decimal integer and LF for each length. Lengths include
the complete patch record through its `END_CHECKPOINT` newline. Readers validate
the file bounds, patch boundaries and final tree marker before publishing state;
editing record widths or newline encoding requires regenerating the index.

MPI runs can set `rebalance_interval_steps` (default 0, disabled) and
`rebalance_minimum_gain` (default 0.1, range [0,1)). On each sampled accepted
step, owner-local reaction, hydro and transport kernel times are summed per
patch, including recursive substeps. Failed retry timings are discarded.
The existing largest-work-first partitioner migrates U/T together only if the
estimated maximum rank load improves by the requested fraction. Geometry stays
replicated; migration uses no interpolation or root gather. Timings exclude
halo exchange, synchronization and migration overhead, so this is not a promise
of wall-clock speedup. Regridding subsequently establishes a fresh cell-weighted
partition; later accepted samples can rebalance it again. These performance
controls are not checkpoint identity: restart obtains fresh timings and may
choose different owners, introducing reduction-order roundoff in budget ledgers.
The serial executable rejects a nonzero rebalance interval.

## Selected mechanisms

When a supported selected bundle is configured, the build also provides
`pelef_reactive_eb_tree_3d_selected` and, with MPI enabled,
`pelef_mpi_reactive_eb_tree_3d_selected`. `selected.nml` is an H2/O2/N2 example
using implicit chemistry. Composition is supplied through `composition_count`,
`composition_species` and `composition_mole_fractions`; fixed-mechanism builds
reject these inputs. The shared source uses the existing mechanism validator,
initial temperature range check and generic source element metadata. Checkpoints
bind bundle SHA-256, composition and integrator, and cumulative chemical
exchange is checked for elemental conservation.

The MPI target links a separate selected runtime without the fixed mechanism
library. It uses the same input and checkpoint identity as the selected serial
target and supports changed-rank and serial/MPI continuation.

The selected executable shares static/scheduled/adaptive controls with the
fixed app. The current focused selected fixture covers a static H2/O2 tree;
it does not qualify detailed-fuel ignition, selected adaptive runs or installation.

## Input

The `eb_tree` namelist uses `root_shape(3)` and `child_count` (0:32).
Children are listed in parent-before-child order: child entry p is tree node
p+1, root is node 1. `parent(p)` names its parent node. `lower(:,p)` and
`upper(:,p)` are inclusive cell indices in that parent; `ratio(p)` defines
refinement. Siblings cannot overlap and each level uses one refinement ratio.
An empty child list runs a root-only tree.

`domain_lower/upper` set the physical domain. `geometry_type`, `plane_axis`,
`plane_position`, `extrusion_axis`, `section_file` and `level_set_file` feed
the existing geometry builder on the full finest lattice. Coarsening must
represent the geometry at every level; unsupported geometry is rejected.

`time_step` is the maximum common interval; `final_time` clips the last step.
With `automatic_timestep=.true.`, `hydro_cfl` (default 0.5) and
`transport_cfl` (default 0.35) bound the common interval using every patch's
cut-cell volume/open-face metrics and global donor signal/diffusivity maxima.
The estimate uses an averaged snapshot, includes hidden parents which are
advanced before restriction, and remains capped by `time_step`/remaining time.
The full-interval diffusion bound is conservatively applied despite split
half steps. It is not a chemistry accuracy estimate.
Rejected steps halve their interval up to `maximum_step_retries` (default 12).
`maximum_steps` bounds accepted steps; reaching it early is an error.
`chemistry_enabled` and `transport_enabled` are independent. Enabled transport
uses viscosity, conduction, species diffusion and barodiffusion.
`chemistry_integrator` accepts explicit/implicit, with `rtol/atol` for reaction.
`riemann_solver` accepts rusanov/hllc/pelec. Initial conditions use the existing
fixed-mechanism density-sheet initializer and its defaults.

## Per-level time stepping

`subcycling=.true.` selects recursive hydro and transport within the existing
R/2 T/2 H T/2 R/2 composition. A hydro child takes `ratio` steps per parent
interval; molecular transport takes `ratio**2`. Ancestor face values are
interpolated at the appropriate level time, and accepted fine transfers are
accumulated and reconciled after each parent interval. Reaction halves still
operate on the exposed composite state around the transport/hydro operators.

The default remains a common interval. `time_step` bounds the root interval;
retries halve it transactionally. In subcycled mode, automatic interval selection
converts each patch bound to the root interval using the cumulative refinement
ratio for hydro and its square for diffusion. Hidden parent updates and global
donor maxima remain included. Checkpoints mark subcycled execution
and reject continuation with a different stepping mode. This composition
retains Euler updates and does not claim second-order time accuracy.

## Scheduled topology changes

`scheduled.nml` demonstrates 5 nodes -> root only -> 4 nodes at times
1.5e-9 and 2.5e-9. `regrid_count` (0:16), increasing positive `regrid_times`,
`regrid_child_count(event)`, `regrid_parent/ratio(:,event)` and
`regrid_lower/upper(:,:,event)` describe complete replacement topologies.
Zero children coarsens to the root. All layouts share root cells and refinement
ratios on common levels. Geometry is built once on the deepest scheduled lattice.

Steps land on event times before conservative remapping and checkpoint saving.
A V3 checkpoint binds the original topology, every scheduled time/topology and
completed stage. Restart reconstructs that stage's geometry before loading
node records, then checks stage/time consistency. Future schedule edits reject
restart. Events beyond the current final time may remain pending.

## Temperature adaptation

`adaptive.nml` evaluates exposed fluid temperatures at the scheduled
`regrid_times`. Set `temperature_adaptive=.true.` and leave prescribed event
boxes/counts at their defaults. `tag_depth` (1:8), `tag_ratios` and
nondecreasing positive `tag_thresholds` define requested refinement levels;
`tag_tile_shape`, `tag_buffer` and `tag_maximum_nodes` (1:33, root included)
control deterministic tiling. The capacity applies before solid-only tiles
are pruned. Exceeding it rejects the adaptation; hot regions are not dropped.

Geometry covers the deepest initial/requested level. Fine hot peaks vote
directly; child-covered parents and solid temperatures do not. Cold states
may coarsen to the root. V4 checkpoints store the complete tagging policy,
evaluation times and all completed generated layouts. Future layouts remain
absent. Restart validates node bounds, parent order, refinement ratios, tile
alignment/extents/order and saved stage/time before continuing. Completed
layouts are restored, without rerunning past temperature decisions.

## Output and restart

`output_file` contains exposed fluid cells only: node, level, local indices,
physical fluid volume, temperature and conserved components. The adjacent
`.budget.csv` gives initial/final composite inventories, boundary, wall,
chemical exchange and residual for each component. Fluid-volume-weighted
cell sums recover the final inventories without counting covered parents.

Set `checkpoint_file`, `checkpoint_interval_steps`, and optionally
`stop_after_checkpoint=.true.`. Continue with `restart_file` and the same
topology, geometry, physics, integrator, time_step, automatic mode and CFL coefficients. Output names and final
limits may change within the existing embedded checkpoint contract. Use a
different checkpoint destination from the restart source.

The application checks the policy envelope, all node records, topology,
geometry, EOS, parent averages, composite budget and EOF before producing
output. Saving serializes to a scratch stream before opening the destination.
Replacement is not filesystem-atomic if an I/O error occurs during copying.

## Scope

This is a serial static, scheduled or temperature-adapted, fixed/selected-mechanism application using common or recursive
Euler operators, spatial PCM donors, hydro outflow and closed diffusion
boundaries. It is not a convergence or physical validation case. MPI and
large-tree scaling remain open for this EB tree application. Plain Euler
retains the small-cell interval restriction; retries do not replace an
accuracy/convergence study.

The EB coefficient default changed to 0.1; old zero-policy checkpoints require
explicit `difmag=0`. See [default-policy checks](../../docs/validation/eb-default-difmag-20260925.md).

### PLM CTU hydro

Set `hydro_method='plm_ctu'` and `riemann_solver='pelec'` for serial or MPI gas
executables, including selected mechanisms. See `plm_ctu.nml`. PCM remains the
default. CTU supports common/recursive intervals and rebuilds its extended
geometry cache after scheduled/adaptive topology changes and restart. Checkpoints
bind the hydro method; changing it on continuation is rejected. The MPI EB spray executable also forwards this selector through its coupled
gas update; its gas checkpoint enforces the same method identity.

`regular_density`, `cut_density` and `initial_pressure` specify the density-sheet
initial state (SI units); defaults remain 0.31, 0.93 and 135000. All must be finite
and positive, and resulting thermodynamic states must satisfy the mechanism EOS.
These values are part of checkpoint identity. Method availability does not imply
that every density discontinuity admits the reconstructed CTU states within the
mechanism temperature bounds; timestep retries may still fail.

## Prescribed inflow and outflow

Serial/MPI fixed and selected gas executables accept explicit `boundary_faces`
in x-low, x-high, y-low, y-high, z-low, z-high order. For example:

```fortran
 boundary_faces='inflow','outflow','outflow','outflow','outflow','outflow'
 inflow_density(1)=.31, inflow_pressure(1)=135000
 inflow_velocity(:,1)=100,0,0
```

Gas executables also accept paired `periodic` labels as described below.
Inlet composition defaults to the initial gas mixture. Override each inflow face
with named mole fractions, for example:

```fortran
 inflow_composition_count(1)=3
 inflow_composition_species(1:3,1)='H2','O2','N2'
 inflow_composition_mole_fractions(1:3,1)=.1,.2,.7
```

Fractions are normalized and converted to mass fractions using the active
mechanism. Names must be unique and present in that mechanism. Count zero
inherits the initial gas; unused entries must be blank/zero. Overrides on
non-inflow faces are rejected. These fields work for fixed and selected
executables. Restart binds the six counts and resolved mass fractions.

`inflow_density` defaults to
0.31 kg/m3, velocity to zero, and pressure must be positive on an inlet.
Leave Mach controls inactive and input `inflow_temperature` at zero; inlet
temperature is recovered from density, pressure and composition. The MPI root
patch owner constructs inlet states once and broadcasts them, including for
idle ranks. Regridding and restart retain those prescribed states.

Both common and subcycled updates forward the controls. Automatic timestep
selection uses boundary donor speeds and diffusivities. Restart context binds
all six labels, densities, pressures and velocity vectors. The MPI spray executable also accepts these explicit boundaries; set
`particle_outflow=.true.` in its spray input. Particle injection remains
controlled by the separate injection schedule. Legacy inputs retain their
existing checkpoint envelope. See the [application checks](../../docs/validation/eb-tree-inflow-application-20260926.md).

## Periodic gas boundaries

For a z-normal plane, x and y can be periodic:

```fortran
 boundary_faces='periodic','periodic','periodic','periodic','outflow','outflow'
```

Both faces of an axis must be periodic and the opposing EB face metrics must
match. Geometry padding wraps periodic donor nodes for PLM/CTU. Common and
recursive gas updates, automatic timestep selection, scheduled regridding and
restart use the same labels. Changing labels on restart is rejected.
MPI EB spray also accepts periodic faces with `particle_outflow=.true.`;
KH child generation and EB deposition/splash with stationary films are supported.
Particle flight, wet-film thickness and film/gas exchange all use the periodic axes.

## Adiabatic domain walls

Fixed/selected serial and MPI gas applications accept `slip_wall` and
`no_slip_wall` in `boundary_faces`. For example:

```fortran
 boundary_faces='slip_wall','slip_wall','no_slip_wall','no_slip_wall','no_slip_wall','no_slip_wall'
```

Wall ghosts reflect momentum. Molecular transport enforces zero mass/species
flux and zero energy flux through stationary adiabatic walls, while retaining
viscous momentum stress. CTU uses reflected configured geometry. Boundary labels
are bound to checkpoint continuation, including across scheduled regridding
and rank-count changes. MPI tree spray accepts these adiabatic wall labels with
`particle_outflow=.true.` and specular particle reflection. KH child generation
is supported at nonperiodic walls. EB impact films can coexist with specular
domain walls; films on the outer domain faces remain unsupported. This input support is not physical wall-flow validation.


## Isothermal domain walls

Fixed/selected serial and MPI gas executables also accept `isothermal_wall`:

```fortran
 boundary_faces=6*'isothermal_wall'
 wall_temperature=6*2000
 transport_enabled=.true.
 automatic_timestep=.true.
```

Temperatures are in kelvin, ordered x-low, x-high, y-low, y-high, z-low, z-high.
Supply six finite positive temperatures when any face is isothermal. Only
isothermal faces use the heat-transfer values; molecular conduction is required.
With no isothermal faces, leave `wall_temperature` at its default zero values.
These stationary walls use no-slip momentum reflection, zero mass/species flux,
and a conductive energy transfer recorded in the boundary budget. Common and
recursive gas splits and automatic thermal timestep selection use these values.
Checkpoint context binds all six temperatures and boundary labels; changing or
removing thermal controls on restart rejects. Nonthermal checkpoint headers
retain their existing format. MPI EB spray also forwards these temperatures to the gas split; use
`particle_outflow=.true.` for explicit boundaries. EB impact films can coexist
with specular domain walls; outer-domain film creation remains unsupported. Nonperiodic thermal walls permit KH child generation. Particles reflect specularly at thermal walls. This is application integration, not physical validation.


Tree checkpoint records support PCM with explicit boundary labels as well as
PLM/CTU. The tree serializer supplies an embedded-patch context while retaining
boundary/solver identity and thermal/inlet validation. This does not enable PCM
physical-boundary stepping in the separate single-level application.
