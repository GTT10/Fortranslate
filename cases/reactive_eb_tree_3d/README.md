# Multilevel EB reacting tree

Build and run the fixed H2/O2 application:

```sh
cmake --build build/debug --target pelef_reactive_eb_tree_3d
./build/debug/pelef_reactive_eb_tree_3d cases/reactive_eb_tree_3d/input.nml
```

The supplied case has a 4-cubed root, two children, and grandchildren with
separate parents (three levels, ratios two). A planar wall cuts z=0.31.
Two 1e-9 steps run reaction, molecular transport and HLLC hydro.

## Reacting mixture and restart acceptance case

`reacting_restart.nml` uses the selected `h2o2_full` mechanism with H2/O2/N2
mole proportions 2:1:3.76. Native MOL with two SDC iterations, molecular
transport and retained internal energy advances a 4-cubed root and one
partial child through three 1e-7 s steps. The child is removed at the first
step and recreated at the second. This is a short reacting integration case,
not an ignition-delay, fixed-PeleC, or production-mesh validation.

With the test suite and MPI enabled, build the selected full-mechanism
fixtures and run the common/recursive acceptance checks:

```sh
cmake --build build/debug --target test_selected_full_reactive_eb_tree_3d test_selected_full_mpi_reactive_eb_tree_3d
ctest --test-dir build/debug -R '^regression_selected_eb_tree_reacting_restart_' --output-on-failure
```

The checker requires measurable hydrogen consumption, a reaction-off
control, nonnegative species, independent H/O/N chemistry inventories,
CSV-integrated gas budgets and mixture NASA7 internal energy. It also checks
serial/2/3-rank final fields, 2-to-3-rank and serial/MPI continuation through
the second regrid, cumulative restart budgets and rejection of changed
chemistry settings. Each invocation writes to a unique directory with input
and binary hashes, commands, logs and `checks.json`; failed attempts remain
distinct from passed acceptance runs.

`reacting_stiff.nml` selects the existing CVODE backend at 250000 Pa with
the same densities, mixture, mesh, transport and 300 ns duration. It needs
`PELEF_ENABLE_SUNDIALS=ON` and the pinned double-precision Fortran SUNDIALS
installation. The initial regular/cut temperatures are approximately
2028/1796 K. The native first-cell diagnostic needed 391 internal steps for
one 100 ns interval; this case exercises the intended stiff-reactor path
without relaxing reaction tolerances.

The focused `regression_selected_eb_tree_hot_cvode_restart` check runs only
the serial case, reaction-off control, two-rank case, checkpoint stop and
three-rank continuation. It reuses the independent inventory/EOS checks and
requires hydrogen consumption above `1e-5` of its initial inventory. The
broader serial-interchange and metadata refusal matrix stays in the existing
tests. The hot check also audits intermediate snapshots and two invalid
output policies without additional time integration.

### Full-coordinate forced CVODE

The opt-in `sdc_reactor='cvode_full'` selects the full species-density and
temperature reactor in fixed/selected serial and MPI EB-tree builds with
SUNDIALS enabled. It uses BDF order at most two, analytic species-density
Jacobian columns and a centered temperature column. It requires reacting SDC;
`native` remains the default and `cvode` retains its reduced-coordinate solver.
The caller supplies tolerances. Endpoint mass admission and transactional
rollback are independent of the requested integration accuracy.

The tree checkpoint stores a distinct reactor identity. Restart with `native`
or `cvode` cannot consume a `cvode_full` checkpoint. The focused
`regression_selected_eb_tree_full_cvode_restart` exercises reacting serial/MPI
execution, changed-rank continuation and reactor mismatch rejection. This
application lifecycle check does not establish fixed-PeleC trajectory or
physical validation; the full-coordinate reactor's earlier 32-cell comparison
belongs to an isolated diagnostic build.

## Fresh MOL reaction history

Fresh reacting gas MOL runs initialize the predictor's reaction history using
private reaction trials at the first attempted timestep. Recursive levels use
their divided intervals. Solid EB cells are excluded, and gas, temperature and
held internal energy are not advanced by initialization. A failed trial leaves
the history unchanged; a rejected first gas step restores it before retrying
with the shortened interval. Restarts retain the checkpoint history, and the
normal accepted-step history takes over after the first step. This change is
connected to gas tree applications; the separate particle startup is unchanged.

## Intermediate gas output

Set `output_interval_steps=10` in `&eb_tree` to save the initial/restored
state, every tenth accepted coarse step, and the final/checkpoint-stop state.
The default, zero, keeps final-only output. The serial and MPI gas
executables share this option; particle/film snapshots are not included and
the spray executable rejects a nonzero interval.

For `output_file='result.csv'`, accepted step 10 writes:

- `result.csv.step10.csv`: exposed composite cells with the same gas,
  species, temperature and available heat-release/UEINT columns as the final CSV.
- `result.csv.step10.csv.budget.csv`: cumulative conservation budgets.
- `result.csv.step10.csv.time.csv`: `step,time_s,nodes,regrid_stage`.

Snapshots are taken after synchronization, rebalancing and regridding, never
from rejected trials. A continuation starts at the checkpoint's absolute
step and physical time; the output interval may change on restart. Filenames
are numbered by step, so sort numerically or use the metadata time. Select a
new output directory/prefix to retain independent runs; matching snapshot
paths are replaced, like the final output. The final `result.csv` and its
budget retain their existing format. This output enables trajectory
comparison but does not itself establish agreement with PeleC.

## Constant gas transport

Serial/MPI gas and MPI spray SDC accept a separate
`&reactive_constant_transport` group in the tree input. Set `enabled=.true.`,
`viscosity` (Pa s), `conductivity` (W/(m K)), optional `bulk_viscosity`
(Pa s), and optional `diffusivity` (rhoD in kg/(m s)). The default
`diffusivity=-1` retains molecular species diffusion. Constants reach common
and recursive stages, physical/embedded walls and automatic timestep bounds.
Restart checks the enabled state and all resolved coefficients. The group
requires enabled SDC transport; HIT scale resolution uses the regular HIT
application. Nonzero Soret controls are rejected by the pinned PeleC contract.

SDC MOL and CTU spatial stages with `none`, `flux_redist` or `state_redist` zero every covered-cell
source and set the density source to the ordered sum of species mass sources,
following pinned `EB.cpp::pc_post_eb_redistribution`. Species, momentum, energy
and the already assembled redistribution registers retain their values.
The application checkpoint tag `PELEF_EB_TREE_SPATIAL_DENSITY_RATE_V2` binds this
spatial-source ordering for every SDC path.

## Method of lines (MOL)

Fixed and selected serial/MPI gas executables accept these additions to the
SDC input:

```fortran
 time_integrator='sdc', mol_mode=.true.,
 hydro_method='plm_ctu', riemann_solver='pelec',
 boundary_faces=6*'outflow',
```

`mol_lax_friedrichs=.true.` selects the LF interface flux; the default MOL
interface is acoustic. `subcycling=.true.` enables recursive level intervals.
MOL and LF choices must match the saved checkpoint. The existing CTU input
selector remains required by the SDC input contract; `mol_mode` selects the
actual hydro update. Three-level serial/MPI comparisons and serial-to-MPI
restart are recorded in [the connection evidence](../../docs/validation/serial-mol-tree-20260927.md).
The nested LF negative-fuel case remains unresolved; these focused checks do
not establish positivity for arbitrary cut-cell states.

MOL peer/coarse-fine synchronization uses the pinned stage registration
weights in both common and recursive intervals. When chemistry is disabled and
`sdc_iterations>1`, the reference registers only the old-time half flux even
though the gas update and physical boundary ledger accept both trapezoidal
halves. The implementation retains that distinction. Fine StateRedist transfers
enter with registered flux; coarse transfers accumulate on every evaluation.
Application checkpoints identify this ordering with
`PELEF_EB_TREE_MOL_REGISTER_REFLUX_V2`; preceding MOL V1 checkpoints are rejected.

For a rotating gas frame, set `rotating_frame_enabled=.true.`,
`rotation_axis=2`, `rotation_angular_speed=120.0`, and
`rotation_axis_point=0.25,0.4,0.5` in `&eb_tree`. The axis uses PeleC's
zero-based x/y/z numbering (0/1/2); speed is in rad/s and the axis point uses
domain coordinates. This requires `mol_mode=.true.`. The application converts
the initial gas energy to the relative-energy convention; inlet states remain
physical and are converted when applied. The frame settings must match the
checkpoint on restart, including a serial/MPI rank-count change. This option
is available in the serial and MPI gas executables and in MPI EB spray SDC,
including recursive particle/film exchange and regular-cell LES. The spray
inventory, injection and particle-outflow ledgers use the rotating-energy convention. The MPI
spray application has a bounded two-rank, two-step check; this does not establish
particle-trajectory parity or physical spray validation.

## Gas LES

With `PELEF_ENABLE_SPRAY=ON`, fixed and selected serial/MPI gas executables
accept a separate namelist in the tree input, without a particle input:

```fortran
&eb_tree_les
 les_enabled=.true.,les_model='smagorinsky',les_cs=.17,
 les_prt=.7,les_sct=.7,les_ci=0,les_cfl=.1
/
```

This requires SDC and completely regular geometry. The algebraic model also
accepts `les_model='wale'` or `'vreman'`. `les_model='dynamic'`
enables `les_test_model`, `les_test_ratio`, `les_grid_ratio` and
`les_harmonic_mean`; the defaults are `'box'`, 2, 1 and `.false.`.
Common and recursive steps use the spatial callback with StateRedist;
Dynamic LES requests the filter's actual donor width. LES also constrains
the timestep. Enable, model, filter and coefficient settings must match a
checkpoint when restarting across rank counts or between serial and MPI.
The runtime remains experimental.

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
`initial_timestep_factor` (default 1) scales the first automatic estimate.
`maximum_timestep_growth` (default 1.1) limits later estimates relative to the
last accepted root interval; the first interval after regridding cannot grow.
Event and final-time caps follow these limits. Automatic mode snaps to the
final time within 0.1% of an interval, matching the reference stop-time rule.
Checkpoints retain both controls,
the last accepted interval and the post-regrid flag; older automatic-mode
checkpoints without that history are rejected.
The full-interval diffusion bound is conservatively applied despite split
half steps. It is not a chemistry accuracy estimate.
Rejected steps halve their interval up to `maximum_step_retries` (default 12).
`maximum_steps` bounds accepted steps; reaching it early is an error.
`chemistry_enabled` and `transport_enabled` are independent. Within enabled
transport, `viscosity_enabled`, `thermal_conduction_enabled`,
`species_diffusion_enabled` and `barodiffusion_enabled` select the physical
operators; all default to true. Barodiffusion is active only with species
diffusion. These controls apply to serial/MPI gas and MPI spray, split/SDC,
common/recursive stepping and automatic timestep estimates. Restart checks
bind the effective settings. Isothermal walls require thermal conduction.
Viscosity includes its energy work and species diffusion includes enthalpy
transport. These are physical-operator controls, whereas pinned PeleC's
`diffuse_*` switches mask output flux components after calculating diffusion;
partial-disable comparisons need to account for that distinction.
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


## Embedded-wall heat and no-slip

Fixed/selected serial and MPI gas SDC applications accept:

```fortran
 time_integrator='sdc', transport_enabled=.true.,
 thermal_conduction_enabled=.true., embedded_wall_temperature=400,
 viscosity_enabled=.true., embedded_noslip=.true.,
```

These options apply to the cut surface inside the domain. Positive temperature
is in kelvin; zero disables prescribed embedded-wall heat transfer. No-slip is
optional and requires viscosity. CTU/MOL and common/recursive SDC account for
wall heat and impulse. The stationary viscous wall performs no mechanical work.
A wall transport frequency caps the time interval even with
`automatic_timestep=.false.`. Both wall settings must match the checkpoint.
See [serial/MPI integration evidence](../../docs/validation/serial-eb-wall-20260927.md)
for the bounded comparisons and restart checks.

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

## SDC common and recursive steps

`redistribution='state_redist'` selects weighted state redistribution for serial/MPI
common and recursive SDC. Its target is `state_redist_target_volume_fraction` (default
0.5, allowed range greater than zero through one). `state_redist_max_order` selects
0 (default), 2 or 4. Orders 2/4 reconstruct at overlap-weighted neighborhood
centroids and preserve component inventories; regular stencils use the selected
limited slope order. `sdc_state_redist.nml` provides a reacting fixture. The method,
target and reconstruction order must match on restart, including rank-changing
continuation. These interfaces do not establish fixed-reference numerical parity. Add `subcycling=.true.` for recursive intervals. See the
[serial coupling/restart checks](../../docs/validation/sdc-eb-tree-state-redistribution-20260927.md)
[MPI common-step checks](../../docs/validation/sdc-mpi-eb-state-redistribution-20260927.md),
and [recursive checks](../../docs/validation/sdc-state-redistribution-recursive-20260927.md).

`redistribution='flux_redist'` enables composite small-cell redistribution for
serial/MPI common and recursive SDC. `sdc_flux_redist.nml` provides a reacting
fixture; add `subcycling=.true.` for recursive intervals. The default is `none`.
Checkpoint continuation requires the same redistribution setting. See the
[serial application/restart checks](../../docs/validation/sdc-eb-tree-redistribution-app-20260927.md)
and [MPI checks](../../docs/validation/sdc-mpi-eb-redistribution-20260927.md).
Recursive continuation is covered by the
[serial checks](../../docs/validation/sdc-eb-redistribution-recursive-20260927.md) and
[MPI checks](../../docs/validation/sdc-mpi-eb-redistribution-recursive-20260927.md).

`sdc.nml` selects `time_integrator='sdc'`, `sdc_iterations=2` and
`external_forcing(3)` (force per fluid volume). This path requires PeleC CTU,
explicit boundary faces. Set `subcycling=.true.` for recursive SDC intervals
with one refinement-ratio worth of substeps per finer level. Scheduled/adaptive regridding also
transfers signed reaction estimates by fluid-volume restriction and PCM
prolongation in the serial application. MPI accepts common and recursive
SDC and preserves owner-local reaction estimates and forcing inventories in V4
checkpoints. Scheduled/adaptive regridding transfers signed
reaction estimates by fluid-volume restriction and PCM prolongation. Load rebalancing
transfers reaction estimates with gas and updates the checkpoint source layout. MPI EB spray SDC remains unavailable.
See the [MPI recursive/restart checks](../../docs/validation/sdc-mpi-eb-subcycle-20260927.md).
See the [serial recursive restart checks](../../docs/validation/sdc-eb-tree-subcycle-20260927.md).
See the [regrid/restart checks](../../docs/validation/sdc-eb-tree-regrid-20260927.md).

Serial fixed/selected source routing preserves reaction estimates and forcing
inventories in V3 tree checkpoint records. The SDC budget adds a `forcing`
column. The fixed executable's continuous/restart comparison is recorded in
[these checks](../../docs/validation/sdc-eb-tree-app-20260927.md).

The same `sdc.nml` fixture runs with `pelef_mpi_reactive_eb_tree_3d`.
A two-rank run and one-rank checkpoint continuation are recorded in the
[MPI application checks](../../docs/validation/sdc-mpi-eb-app-20260927.md).

## Optional momentum limiter for coarse halos

`separate_momentum=.true.` selects an experimental coarse-halo interpolation
policy: each momentum component has its own bounded slope, while density,
species and energy retain a shared scalar limiter. The default is false.
It requires `time_integrator='sdc'`, `mol_mode=.true.`,
`dual_energy_enabled=.true.` and `redistribution='none'`. Common and recursive
stepping accept it; spray and enabled LES are excluded.

MPI ranks must agree on the choice. A checkpoint written with this option
requires it on restart; switching either direction is rejected before state
loading. The default checkpoint envelope is unchanged. Kernel/operator and
native source/history replay tests pass. Selected reacting common/recursive
application lifecycle checks also pass; see the
[application evidence](../../docs/validation/eb-momentum-application-20261021.md). This option changes the numerical policy relative to the original
pinned reference and is not physical validation. See the
[design and evidence](../../docs/design_decisions/eb-halo-momentum-limiter.md).

### Density-coupled velocity option

`consistent_velocity=.true.` together with `separate_momentum=.true.` selects
an alternative momentum slope: the donor velocity multiplies the limited
density slope, and a separately limited relative-momentum slope is added.
This preserves constant velocity when the scalar limiter clips density.
Both options default to false. It uses the same held-MOL, no-redistribution,
no-LES/no-spray admission rules above.

The choice is passed through serial/MPI hydro and transport halos, common and
recursive stages, retry/source-prefix identity and history replay. MPI ranks
must agree; checkpoints record the velocity policy and refuse switching it
on restart. This connects the velocity candidate to ordinary application
stepping. The separate diagnostic pressure-reconstruction candidate is not
selected by this input. See the
[velocity application record](../../docs/validation/eb-velocity-application-20261022.md).

### Paired pressure option

For common or recursive gas held MOL, select all three policies:

```fortran
 separate_momentum=.true., consistent_velocity=.true., paired_pressure=.true.,
 subcycling=.false., redistribution='none',
```

The pressure option uses fluid-volume-restricted fine EOS pressure on both
covered coarse spatial states and incoming fine halos. It requires the
velocity option, held MOL and a nonrotating frame. Set `subcycling=.true.`
for recursive clocks with endpoint-pressure ancestor interpolation.
LES and spray are unsupported and refused. Checkpoints bind the choice;
switching it on restart is refused. It does not change accepted input leaf
states directly or eliminate pressure error from general thermal mixing.
See [the numerical policy](../../docs/design_decisions/eb-paired-pressure.md)
and [recursive application checks](../../docs/validation/eb-pressure-recursive-20261022.md).

### Connected-direction halo option

`connected_directions=.true.` enables interpolation along connected coarse
stencil directions beside a covered neighbor. It defaults to false and
requires `separate_momentum=.true.` and its held MOL/SDC admission rules
(no redistribution, LES or spray). It may be combined with
`consistent_velocity` and `paired_pressure`.

Blocked directions contribute zero slope. Disconnected diagonals in the
retained stencil force PCM. Serial/MPI common and recursive dispatch pass
the option, and MPI ranks must agree. Checkpoints bind it with
`PELEF_EB_TREE_CONNECTED_DIRECTIONS_V1`; switching it on restart is refused.
The integration verification is recorded separately from analytic-shear
convergence in [the application record](../../docs/validation/eb-connected-application-20261022.md).
Local transverse-error convergence and general geometry flow accuracy remain open.

## SDC CVODE reactor

The gas EB-tree `&eb_tree` namelist also accepts
`pinned_eos_convergence=.true.` when `pinned_eos_extrapolation=.true.`.
This opt-in caloric inverse uses a 1e-12 K temperature correction threshold
and up to 5000 iterations; the default remains false (1e-6 K, 200 iterations).
It does not change chemistry tolerances. Serial/MPI fixed and selected
frontends pass the policy to the existing thermodynamic routines.
Checkpoint policy identity includes this setting, so continuation requires
the same value. MPI rejects differing policies across ranks before gas
initialization. The selected full-mechanism gas reacting/regrid/restart check
passes; see the [application record](../../docs/validation/eb-eos-application-20261021.md);
this is not spray or LES qualification.

SUNDIALS-enabled fixed/selected serial and MPI gas executables, and MPI particle
executables, accept `sdc_reactor='cvode'` in `&eb_tree`. It requires reacting
SDC. `native` remains the default; this selector is independent of the split-step
`chemistry_integrator`. Common and recursive stepping forward the selected
reactor, including the particle path's optional LES/film callbacks.

Restart requires the same reactor. Native checkpoint control values remain
unchanged; CVODE uses an additional control flag that older readers reject.
See [focused gas/particle checks](../../docs/validation/eb-tree-particle-cvode-20260927.md)
for the tested combinations and limits.

The opt-in `dual_energy_enabled=.true.` SDC gas path retains a separate
`rho*e` (UEINT) field through timestep selection, regrid, owner rebalance and
serial/MPI checkpoint continuation. It connects MOL or Godunov CTU on common
or recursive steps with no redistribution or prepared StateRedist. Reset
controls are `allow_small_energy`, `allow_negative_energy`,
`update_total_from_internal` (all default true) and `dual_energy_eta2` (default
`1e-4`); checkpoints bind all four controls.
CSV output adds `rho_internal_energy_J_m3`. The budget separates
`numerical_energy` from chemistry, forcing and physical boundary work.
`internal_energy_pressure_work` reports cumulative signed Godunov `-p div(u)`
volume work on the energy row; it is an UEINT diagnostic and is not subtracted
from the UEDEN conservation residual. MOL reports zero for this Godunov term.
Held checkpoints use serial V9 or indexed MPI V10 to retain this scalar history
across rank changes. The readers retain V7/V8 compatibility for MOL; Godunov
restart requires the pressure-work history in the new format.
MOL LES and embedded viscous/thermal wall bindings have focused checks.
Godunov reacting/reset/LES trajectories and paired reactive parent replay need
broader qualification. Held particle/film coupling remains unfinished and is
rejected before initialization. All remain part of the port's completion scope.
Existing gas-only input/output formats are retained when the option is off.

## Explicit gas initial state

`initial_state_file` optionally loads `PELEF_EB_TREE_STATE_V1` for a fresh,
nonrotating gas run with `dual_energy_enabled=.true.`. It conflicts with
`restart_file`; the default remains density-sheet initialization. MPI rank 0
reads and broadcasts the records. Candidate owner-local arrays are committed
only after the complete input passes admission.

The text layout is: magic line; node and species counts; one species name per
line in mechanism order; then, for each node, its one-based index and three
shape dimensions followed by all cells in k/j/i loop order (i fastest).
Each cell row contains `node i j k`, the conserved vector
`rho rhou rhov rhow rhoE rhoY(1:nspecies)`, held internal-energy density,
and temperature. Covered coarse donors and solid cells must be included.
Nonfinite values, negative species, inconsistent mass/energy/temperature,
wrong order, missing cells and trailing records are refused. Empty cells
are allowed only where geometry has no fluid. File values are SI.
For uncovered cells, total minus kinetic energy must equal held internal
energy. Covered coarse donors may retain nonnegative unresolved kinetic
energy from independent conservative restriction of total and held energy;
a negative deficit is refused. Their held energy must still admit EOS
inversion; a positive independently restricted temperature is retained as
the inversion seed. Uncovered cells require temperature agreement with the
EOS. The importer preserves these donor
values; it does not assert that supplied coarse values are fine-cell averages.

This input path does not itself establish fixed-PeleC agreement. The geometry,
transport, timestep, reaction tolerances and numerical policies still come
from the application configuration and must match the reference comparison.

### Ordinary MOL transport thermodynamics

Ordinary nonrotating MOL uses total energy minus resolved kinetic energy to
recover temperature and pressure for normal-face transport, consistently with
its hydro fluxes and the pinned PeleC path. After restriction, this can differ
from the separately retained internal energy. The retained scalar continues to
follow its own evolution; general transport, rotating frames, experimental
separate-momentum halo policies and wall-source recovery retain their existing
contracts. Serial and MPI use the same selection. The dated
[transport recovery evidence](../../docs/validation/eb-full-cvode-application-20261024.md)
separates bounded reference comparisons, restart equivalence and open C2 work.

### Selected-cell MOL diagnostics

Serial held-energy SDC MOL with subcycling and redistribution `none` accepts
`mol_trace_count` (0 by
default, maximum 16), `mol_trace_nodes`, and `mol_trace_cells(:,probe)`.
For example, `mol_trace_count=1, mol_trace_nodes(1)=2,
mol_trace_cells(:,1)=3,3,8` records one fine cell. MPI requests are rejected
explicitly; this diagnostic is not an MPI trace implementation.

The application writes `<output_file>.mol_stage.csv` and
`<output_file>.mol_sync.csv`. Component indices follow conserved-state order:
density, three momenta, total energy, then species densities in mechanism
order. Stage tables contain cell and reconstructed face pressures, raw and
centroid hydro fluxes, supplied tree flux divided by the interval, hydro/wall
rates and the observed state. Sync tables retain interval-integrated fluxes,
state, and raw/applied closure increments. Each has node/cell/face, depth,
interval-start time, interval, step and event identifiers. These are read-only
attempt diagnostics, including any subsequently rejected attempts; they are
not accepted-step receipts. Hydro rates exclude reaction/transport, so their
balance alone does not reconstruct the full reacting SDC trajectory. Consult
face aperture and volume fraction, especially for closed faces and cut cells.

Reacting runs additionally write `<output_file>.mol_reaction.csv` after each
successful forced-reactor pass, before the subsequent energy recovery/reset.
This table records the interval's old state, updated state, actual forcing
rate, old/new spatial and external rates, and the returned reaction source.
The forcing includes the hydro contribution and any applicable predictor
correction; it must not be replaced by the hydro-only stage flux divergence.
For species, the reaction source is the reactor endpoint minus its forced
nonreacting endpoint, divided by the interval. Its energy entry is the
returned total-energy projection correction, not chemical heat release.
Density and momentum entries are zero; density can retain roundoff from
reconstructing the sum of species. As with the other tables, rows describe
attempts and are not proof that a step was accepted.

Stage rows also include `rho_left/right`, `normal_velocity_left/right`, and
`sound_left/right` in SI units. Density and velocity are the successful
Riemann call's reconstructed states; normal velocity follows the face axis.
Sound speeds are the two donor cell-center values whose mean is supplied to
the PeleC Riemann path. Closed faces retain zero diagnostics and must be
excluded using aperture. These fields do not include transverse velocity or
species composition and therefore are not a complete arbitrary-state replay.

The opt-in trace also writes `<output_file>.mol_stencil.csv`: four prepared
primitive donor cells per face, before characteristic reconstruction. Offsets
-2,-1,0,1 are relative to the right donor; -1 and 0 are the left/right cells.
Columns contain SI density, three velocities, pressure, donor sound speed,
and `Y_1` through `Y_n` in mechanism order. Covered cells retain zero primitive
values. This table includes halo inputs but does not record EB connectivity;
use the face geometry when interpreting or replaying it. Allocation and output
occur only when MOL diagnostics are requested.
