# Experimental regular 3D spray combustion (0.246.0)

This is a working opt-in Eulerian--Lagrangian path, not production diesel-spray
qualification. It couples the existing regular 3D gas solver to
single-component point parcels, molecular transport, a deviatoric Smagorinsky
closure and native or CVODE cell chemistry. It does not replace the existing
AMR/EB/MPI applications or claim those features are connected to particles.

## Build and run

A fixed H2/O2 gas mechanism with water droplets requires no new dependencies:

```sh
cmake -S . -B build/spray -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DPELEF_ENABLE_TESTS=ON -DPELEF_ENABLE_SPRAY=ON
cmake --build build/spray --parallel 4
ctest --test-dir build/spray -L spray --output-on-failure --no-tests=error
mkdir water-run && cd water-run
../build/spray/pelef_spray_3d ../cases/spray_3d/water.nml
```

For liquid-fuel combustion, the reference fixture is FFCM1_Red (21 species,
124 reactions) from the **already pinned** PelePhysics revision. The helper
checks the upstream YAML and license bytes before importing. It requires
Cantera 3.2.0 and never overwrites different local source data.

```sh
python3 -m pip install cantera==3.2.0
python3 tools/prepare_spray_reference.py --output-dir build/methanol-reference
cmake -S . -B build/spray-cvode -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DPELEF_ENABLE_TESTS=ON -DPELEF_ENABLE_SPRAY=ON \
  -DPELEF_ENABLE_SUNDIALS=ON -DPELEF_ENABLE_CANTERA_REFERENCE=ON \
  -DCMAKE_PREFIX_PATH=/absolute/path/to/your/sundials/install \
  -DPELEF_MECHANISM_BUNDLE="$PWD/build/methanol-reference/ffcm1_red.json" \
  -DPELEF_MECHANISM_SOURCE="$PWD/build/methanol-reference/ffcm1_red.yaml"
cmake --build build/spray-cvode --parallel 4
ctest --test-dir build/spray-cvode -L spray --output-on-failure --no-tests=error
mkdir methanol-run && cd methanol-run
../build/spray-cvode/pelef_spray_3d_selected ../cases/spray_3d/methanol_burn.nml
```

SUNDIALS must provide the static double-precision official Fortran CVODE,
serial-vector, dense-matrix and dense-linear-solver interfaces built with a
compatible Fortran compiler. The qualification workflow pins SUNDIALS 7.2.0 at
`71a4cc9ad5e7bc8b4e33a1ca9795b4e96883f9a6`. `chemistry_method='cvode'` fails
explicitly when that backend was not built. With spray disabled (the default),
the preceding application/install target inventory remains unchanged.

Offline reference preparation uses both `--yaml PATH --license PATH` with the
same exact pinned checks. Source/license hashes and the illustrative liquid
parameters are recorded in [the manifest](../references/spray_methanol.json).
The reference helper retains the upstream license; it does not select a license
for this project.

## What the supplied cases establish

| Input or test | Purpose | Not established |
| --- | --- | --- |
| `water.nml` | Water evaporation and two-way gas feedback without reaction | Diesel fuel properties or injector accuracy |
| `methanol_burn.nml` | Closed uniform liquid-fuel evaporation through oxidation and net heat release | A spatial spray flame, lift-off, or measured ignition delay |
| `methanol_cone.nml` | Nonuniform cone, scheduled injection, transport, reacting cells and LES together | A physical inflow/nozzle boundary or turbulence statistics |
| `check_spray_cantera.py` gas-only cases | Independent full temperature history and final composition against Cantera | An external reference for the liquid model |
| `check_spray_cone.py` | All-phase conservation, elements and exact coupled restart with injection/chemistry/LES | Experimental spray validation or parallel scalability |

For a new run, output paths and checkpoint destinations must not exist. Outputs
are `.history.csv`, `.gas.csv` and `.parcels.csv`. History records gas/liquid
mass, system momentum and total energy, injection mass, balance residuals,
temperature bounds, live parcels and Sauter mean diameter. A parcel's `mass`
is the mass of **one** droplet; its multiplicity multiplies all transfer terms.

## Models and conventions

SI units are used throughout. Each liquid population is dilute, spherical, single-component,
spatially isothermal, with constant or temperature-dependent polynomial density. `vapor_name` must identify a gas
species in the selected mechanism. Liquid property numbers in the examples
are illustrative constants, **not** a fitted methanol/water property database.

The reference energy is `e_l(Tref) = h_v(Tref) - Lref`, using the selected
mechanism's NASA7 gas enthalpy. Thereafter `e_l(T) = e_l(Tref)+integral(cp_l(T), Tref, T)` and
the latent heat in the droplet energy equation is `h_v(T)-e_l(T)`. The
incompressible liquid pressure-work correction and surface energy are neglected.
This choice aligns chemical formation-energy references across phases. The
gas receives minus the parcel's discrete change of mass, momentum and **total**
energy; latent heat, evaporated kinetic energy and drag work are not added a
second time.

The surface saturation pressure uses a **constant-Lref** Clausius--Clapeyron
law anchored at `(Tref, pref)` by default. States with `psat >= pgas` are rejected
by default; the optional surface cap below limits the vapor-state evaluation temperature. Condensation is not modeled;
a saturated or supersaturated ambient produces zero evaporation. Film state
uses the one-third rule, NASA7 thermodynamics and existing molecular transport.

The vapor/carrier diffusivity is the carrier-normalized harmonic binary value,
`D_vc = sum(j!=v, X_j) / sum(j!=v, X_j/D_vj)`. It is evaluated directly from the
binary matrix and carrier composition, recovering the exact binary coefficient
for a single carrier. Do **not** apply a PelePhysics `rho*D` conversion to
PeleF's differently defined mixture-averaged coefficient. A regression reproduces
the erroneous molecular-weight conversion and guards the corrected binary limit.

Drag uses the [Schiller--Naumann sphere correlation](https://cpp.openfoam.org/dev/SchillerNaumann_8C_source.html) (`Cd*Re=24*(1+0.15*Re^0.687)`
below Re=1000; `Cd=0.44` above). Transfer supports Ranz--Marshall or an
Abramzon--Sirignano-style blowing correction. The latter iterates Nu/Bt,
uses the film mixture heat capacity, and fails on nonconvergence. See the
[PeleMP equations](https://amrex-combustion.github.io/PeleMP/Equations.html) and
[Abramzon and Sirignano (1989)](https://doi.org/10.1016/0017-9310(89)90043-4).
The constant-latent saturation law, drag correlation and regular-grid CIC
source distribution differ from PeleMP: this is not bitwise PeleMP parity.

Frozen-film drag and heat relaxation are integrated exponentially, evaporation
uses a nonnegative diameter-squared update, and parcel steps are limited by
relative changes, at most a small temperature change, and cell traversal.
A minimum-diameter event transfers the **entire** residual mass/momentum/energy.
By default particle locations wrap periodically. Gas interpolation and source deposition
use the same eight cloud-in-cell weights. Inadmissible gas or liquid candidates
are rejected with step reduction; exhausting the retry/substep limit rolls back
the entire caller-visible coupled step. Sequential parcel feedback is
deterministic but not parcel-order-independent or MPI distributed.

The SGS model uses `mu_t=rho*(Cs*Delta)^2*sqrt(2*S:S)`, deviatoric stress,
turbulent heat conduction and corrected species diffusion with zero net species
mass flux. Conservative face energy flux includes stress work and species
enthalpy diffusion. Its SSPRK2 update checks positivity/EOS and a diffusive
step bound. `Cs=0` is an exact no-op. There is no isotropic SGS kinetic-energy
model, wall damping, dynamic coefficient, stochastic parcel dispersion or
unresolved turbulence--chemistry interaction closure. See
[the PeleC LES discussion](https://amrex-combustion.github.io/PeleC/LES.html) for
context; this smaller closure is independently implemented and is not all of
PeleC's compressible LES model.

A coupled step applies half spray/LES/chemistry updates, the existing gas
hydro/molecular step, and reverse half updates. Each CVODE cell call uses that
cell's current density/internal energy and releases its context even on failure.
It does not reuse an integrator with stale density after evaporation. All gas,
temperature and parcel changes are committed together. Symmetric operator
ordering does **not** make the frozen-film parcel scheme second-order; temporal
and mesh convergence must be checked for each research problem.
`check_spray_convergence.py` checks history-error reduction at three outer
timesteps for the closed methanol example; it does not certify all meshes or
second-order spray accuracy.

## Injection and restart

The cone uses deterministic equal-solid-angle sampling with an azimuth based on
the monotonically assigned parcel ID. Droplets share the input temperature. `diameter_shape=0` (default) is
monodisperse at `drop_diameter`. A positive shape in [0.25,20] enables a
truncated Weibull/Rosin--Rammler **mass CDF** between `diameter_min` and
`diameter_max`, with `drop_diameter` as the scale (not the Sauter mean).
Midpoint inverse-CDF quantiles each carry equal parcel mass; multiplicity
therefore varies as diameter to the minus third power. This reproduces the
specified mass distribution, not a number-based Weibull distribution.
It is not OpenFOAM's differently interpreted `massRosinRammler` correction.
[The underlying truncated distribution](https://api.openfoam.com/2412/RosinRammler_8H_source.html)
is sampled deterministically; it is not a turbulent dispersion model. Prescribed mass flow is introduced in timed pulses;
the final pulse uses the remaining window duration so injected mass is exact.
The fluid time step is clipped to injection events. This is a prescribed
particle source in a periodic box, **not** a resolved nozzle or open inlet.

`stop_after_steps` plus `checkpoint_file` writes a coupled checkpoint;
`restart_file` resumes into a **new** output prefix. The file contains gas and
liquid state, time/step/event counters and original/injected budgets. It embeds
the mechanism identity, application version, backend availability and physical
input context; changing physics or timestep controls is refused. Output paths,
stop controls and a later final time may change. Ordered IDs, shape, finite
values, state/temperature consistency, event counts and balance ledgers are
checked before changing caller state or opening output files.

Publication closes a temporary file before an atomic no-overwrite POSIX hard
link in the destination directory. Existing checkpoints are never replaced.
This requires a Linux/POSIX filesystem supporting hard links; it is not an
`fsync`-durable guarantee against power loss. Files are strict formatted records,
not cryptographically authenticated data. Malformed/truncated/extra/null records
and changed physical context are rejected. Read failure leaves existing state
unchanged. Restart compatibility is intentionally limited to the same recorded
solver version/physics, not an indefinite checkpoint ABI.

## Explicit remaining work

Production use still needs calibrated temperature-dependent liquid properties,
boiling/real-fluid handling, multicomponent liquids, additional breakup/dispersion models and collision, calibrated size distributions
and turbulent dispersion, physical injector and wall-impact models, particle
migration with MPI/AMR/EB, scalable output/checkpointing, detailed diesel-fuel
mechanisms and performance tests. Selected chemistry uses a 2--128-species
limit and currently supported reaction forms; FFCM1_Red does not represent
HXN/HMN/AMN. No claim of measured ignition, penetration or liquid length is made.


## Regular-grid wall and outflow boundaries

`boundary_faces(6)` in `&spray` uses x-low, x-high, y-low, y-high, z-low,
z-high order. The default remains periodic. Supported alternatives are
`inflow`, `outflow`, `slip_wall`, `no_slip_wall`, and `isothermal_wall` for gas; parcels
reflect specularly at each wall type. Opposing periodic faces must be paired.
This ideal reflection is not a sticking, splash, film or impact model.

Physical faces use a clamped, normalized interpolation/deposition stencil with
unique donor cells; particles never sample the opposite side of a physical wall.
Outflow removes a parcel's remaining mass, momentum and energy from the system;
wall reflection transfers normal momentum without changing parcel energy.
The split step accumulates gas surface fluxes and parcel exchange in one ledger.
Physical-boundary history adds `boundary_mass`, `boundary_px`, `boundary_py`,
`boundary_pz`, `boundary_energy`; positive values enter the simulated system.
The reported conservation errors subtract this ledger and injected inventory.

`wall_reflection.nml` supplies an executable example. Change its faces to
`outflow` for particle escape. Checkpoint V2 preserves boundary exchange, and
its context binds all six face types. Periodic runs retain checkpoint V1 and
their existing history columns. Restart reproduces gas/parcel CSVs and the final
history row exactly in the focused [boundary checks](validation/spray-boundaries-20260915.md).

Nonperiodic faces also support the deviatoric Smagorinsky model. Physical-face
stress, heat and species fluxes use mirrored/clamped ghost states; the surface
ledger includes both LES RK stages in both split half-steps. The diffusion bound
uses only faces touching physical cells, excluding the artificial outer ghost
boundary. Set `smagorinsky_constant=.17` and `shear_velocity=10` in the example
to exercise nonzero SGS. See [focused SGS evidence](validation/sgs-boundaries-20260915.md).

This uses the existing constant-coefficient closure without wall damping or a
wall function. Near-wall turbulence qualification, characteristic inlet models,
AMR/EB/MPI parcels and physical wall-impact models
remain open. Injection is still an internal parcel source.


## ETAB representative-drop breakup

Set `breakup_model='etab'`, `surface_tension` (N/m), and `liquid_viscosity`
(Pa s); both material properties must be positive. The default is `none`.
`etab.nml` exercises constant illustrative water properties at high gas density
and relative velocity. These values are a numerical fixture, not a fitted fuel
property model.

The deformation oscillator and ETAB child-size law follow the pinned
[PelePhysics TABBreakup.H](https://github.com/Pele-Suite/PelePhysics/blob/6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b/Source/Spray/BreakupSplash/TABBreakup.H).
Breakup decreases a representative droplet's mass and increases its multiplicity
at constant parcel inventory, velocity and temperature. It preserves system
mass, momentum and the existing thermal/kinetic energy convention, which omits
surface energy. Child radii are bounded by `minimum_diameter/2`; an excessive
substep count or invalid state rejects the complete coupling step.

`deformation` and `deformation_rate` are included in breakup-enabled parcel CSVs
and checkpoint V3 (periodic) or V4 (with physical-boundary inventory). Material
parameters and model identity are bound into checkpoint context. Models without
breakup retain V1/V2. The damped oscillator advances only the remaining interval
at each internal step; zero-slip and overdamped states avoid singular formulas.

The representative-size/deformation update can also partition computational
parcels as described below. Stochastic transverse child velocities, primary atomization,
collision/coalescence and experimental breakup qualification remain open. Source hashes,
license and deliberate scope differences are recorded in
[the reference manifest](../references/spray_etab.json) and
[focused evidence](validation/etab-20260915.md).


## Mixed-liquid thermodynamic API

`spray_mixture_mod` provides the shared thermodynamic operations for up to 16
components of one liquid droplet, distinct from the existing per-port pure liquids.
It accepts pure-component property records and normalized mass fractions.
Density uses additive specific volume, `1/rho = sum(Y_i/rho_i)`; heat capacity
and stored liquid energy use `sum(Y_i*cp_i)` and `sum(Y_i*e_i)`. Each component's
energy retains its mechanism formation reference. No excess mixing energy is added.

`liquid_mixture_properties` evaluates these quantities at a common temperature.
`liquid_mixture_temperature` inverts the energy with bracketed Newton iteration,
inside the active components' common density/heat-capacity fit ranges. Zero-mass
components do not constrain the temperature range. Non-normalized fractions,
duplicate component species, invalid properties or unreachable energies reject
without replacing the outputs.

`advance_mixture_heating` applies a supplied parcel heat amount at fixed composition,
recovers temperature and returns the opposite actual energy increment for the gas.
Failure preserves the parcel and returns zero transfer. The transport update
below supplies mixture heat transfer and component evaporation.
The CLI supports shared or per-port mixed-liquid compositions using the inputs below. The
mixture grid API retains component masses and supports their checkpoint persistence
and combined hydro/chemistry stepping. The
reference's optional boiling-temperature property caps are not applied here.
See [source](../references/spray_mixture.json) and
[focused evidence](validation/spray-mixture-thermo-20260916.md).

### Mixed-liquid surface and component exchange

`liquid_mixture_surface` converts mass fractions to liquid mole fractions and
returns equilibrium surface fractions `X_liquid_i * psat_i / P` (Raoult's law).
It rejects a saturated combined surface and leaves outputs unchanged on failure.
A component's pure saturation pressure may exceed ambient pressure when dilution
keeps its partial pressure below ambient. Optional per-component surface caps
use the actual total gas pressure before applying the liquid mole fraction.
This surface routine returns equilibrium fractions; the film builder below adds
ambient-vapor suppression and transport coefficients.

`advance_mixture_evaporation` takes a parcel, component masses **per physical
droplet**, requested evaporated masses on the same basis, component vapor
enthalpies in J/kg, and heat supplied by the gas to the **whole parcel** in J.
It updates the supplied component-mass array and total droplet mass, recovers the
remaining mixture temperature, and returns gas conserved-inventory increments
(mass, three momenta, total energy, followed by gas species masses).
Components map through their pure-liquid vapor indices; duplicate mappings reject.

For a surviving droplet, the target liquid energy equals old liquid energy plus
supplied heat minus evaporated vapor enthalpy. The gas receives the actual
opposite discrete liquid inventory change, including evaporated kinetic energy;
no separate latent-heat term should be applied again. Velocity and position do
not change in this exchange operation. A completely evaporated droplet returns
its finite entire inventory once, following the existing terminal-drop convention;
there is no residual liquid temperature to solve or additional heat debit.
Overdraw, invalid composition, unsupported energy states and invalid species
indices reject without modifying parcel/component masses and return zero exchange.
This API consumes calculated evaporation amounts; the frozen-gas update below
supplies them. The grid API below retains component masses alongside the parcel array.
See [component-exchange evidence](validation/spray-mixture-exchange-20260916.md).

### Mixed-liquid film transport and frozen-gas update

`build_mixture_film` combines the liquid surface fractions with the ambient gas.
A component with surface mole fraction no greater than its ambient fraction has
zero outward flux. The remaining gas species provide the carrier molecular weight;
the one-third film rule supplies temperature and composition. Each vapor uses the
existing carrier-mole harmonic binary diffusivity. Effective diffusivity is the
surface-mole-weighted average; component evaporation shares are proportional to
`X_surface_i * D_i`, following the reference's component weighting structure.

`mixture_transfer_rates` applies the existing selectable drag and
Ranz–Marshall/Abramzon–Sirignano correlations using the actual mixed-liquid density
and heat capacity. It returns inverse drag time, heat coefficient in W/K,
component mass rates in kg/s **per physical droplet**, and vapor enthalpies.
It uses PeleF's carrier molecular weight, binary diffusivities and film-cp heat
correction; this is not exact upstream coefficient parity.

`advance_mixture_frozen_gas` combines exponential drag, trapezoidal flight,
heat transfer and component evaporation over a supplied interval. Gas state is
held fixed during this operation; returned conserved increments include drag
work and species exchange. It changes the parcel and external component-mass
array only on success. Overdraw or inadmissible temperatures reject the interval;
the caller must shorten it. Evaporating trace components at or below a mass
fraction of 1e-12, or with a predicted residual below that fraction of the current
droplet mass, release their entire residual inventory. They no longer constrain
the grid step through a relative component-mass limit. Suppressed components are
retained. If evaporation reduces the predicted mixture diameter below
`minimum_diameter`, every component is transferred to the gas using the terminal
finite-inventory convention. These transfers conserve mass and energy rather
than discarding small species and renormalizing their survivors.
A zero interval preserves state exactly. This remains
a composable API used by the grid and wall coupling below.
See [rate and update evidence](validation/spray-mixture-rates-20260916.md).

### Mixed-liquid grid coupling and composition restart

`advance_mixture_coupling_3d` advances a parcel array and a persistent
`component_mass(component,parcel)` array together with a regular 3D gas grid.
All parcels use the supplied ordered component table. Each column stores masses
per physical droplet and must sum to that parcel's droplet mass. The routine
samples gas with the existing particle stencil, limits drag/temperature/component
changes and particle travel, and halves rejected feedback intervals. Gas,
particle, component and boundary states commit only after the full update succeeds.

The boundary API supports periodic wrapping, specular wall reflection and
outflow/inflow removal. Boundary budgets include each component's species mass;
`mixture_system_totals` includes both gas and liquid inventories on the same basis.
ETAB/KHRT representative-size updates preserve composition while scaling droplet
mass and multiplicity. Their surface tension and liquid viscosity must be supplied
as mixture properties through the existing options. Computational splitting copies
component masses per droplet into its descendants. Optional KH shed-mass children
inherit the parent composition, with conservative transverse kinetic energy debited
from gas. Children advance through the remaining step; capacity or gas-state failure
rolls back the entire update. KH accumulators and composition survive computational
splitting and restart. Wall-impact/film conversion is available through the optional coupled arguments below.
See [mixed KH evidence](validation/spray-mixture-kh-20260916.md).

Both checkpoint routines accept optional `component_mass`. The `_MIX_V1` suffix
adds a component count after grid dimensions and a component-mass row after each
parcel record. Counts, shapes, nonnegative finite masses and parcel-mass closure
are checked. A reader without composition output rejects this format; invalid
reads preserve existing state. The caller's restart context must include `mixture_evolution_id` and identify the
ordered component properties and mixture settings, as with existing pure-liquid
checkpoint contexts. The reader commits the allocated component matrix together
with the gas and parcel state.

These regular-grid coupling/checkpoint APIs are used by the mixture namelist
and scheduled injection path below.
See [grid/restart evidence](validation/spray-mixture-grid-20260916.md).

### Mixed-liquid reacting-flow step

`advance_spray_reactive_3d` accepts the optional pair `mixture_components` and
`component_mass`. It advances the same mixture state in both spray half-steps,
with the existing hydro, molecular transport, LES and chemistry operators between
them. External acceleration retains its impulse/work ledger. Gas, particles and
component masses commit together only after every stage succeeds; late chemistry
or transport rejection returns zero external/boundary ledgers and preserves inputs.

The mixture pair cannot be combined with the per-port pure-liquid property tables.
Mixed particle wall impact requires `wall_films` and `film_component_mass` together.
Gas wall conditions remain available. No new chemistry
backend or time-order accuracy claim is introduced. The focused check exercises
the native implicit backend, not a new CVODE qualification. See
[reacting-mixture evidence](validation/spray-mixture-reactive-20260916.md).

### Mixed-liquid application inputs

Set `mixture_count` to 1–16 and supply that many entries in:

| Input | Meaning / units |
| --- | --- |
| `mixture_vapor_names` | Distinct gas species names |
| `mixture_fractions` | Initial liquid mass fractions, nonnegative and summing to one |
| `mixture_densities` | Constant-model component densities, kg/m³ |
| `mixture_cps` | Constant-model component heat capacities, J/(kg K) |
| `mixture_reference_temperatures` | Component reference temperatures, K |
| `mixture_reference_pressures` | Reference saturation pressures, Pa; default 101325 |
| `mixture_latent_heats` | Reference latent heats, J/kg |

These inputs define the component formation-referenced energy and Clausius vapor
pressure. Unused entries must remain at their defaults. `mixture_count=0` keeps
the existing pure-liquid mode. Component property models can be set independently:

| Input | Values / layout |
| --- | --- |
| `mixture_density_models`, `mixture_cp_models` | `constant` (default) or `polynomial` |
| `mixture_density_coefficients`, `mixture_cp_coefficients` | `(4,16)` cubic coefficients, `a + T*(b + T*(c + T*d))`, SI units, T in K |
| `mixture_density_temperature_ranges`, `mixture_cp_temperature_ranges` | `(2,16)` valid temperature bounds in K |
| `mixture_vapor_pressure_models` | `clausius` (default) or `antoine` |
| `mixture_antoine_coefficients` | `(4,16)` values A,B,C,D for `psat = D*10**(A-B/(T+C))`, Pa and K |
| `mixture_antoine_temperature_ranges` | `(2,16)` valid temperature bounds in K |
| `mixture_boiling_models` | `reject` (default) or `surface-cap`, the latter requires Clausius |
| `mixture_critical_temperatures` | Component critical temperatures in K for surface-cap |

For example, `mixture_cp_coefficients(:,1)=4000,1,0,0` selects the first
component's coefficient column. Coefficients and ranges for inactive models must
remain zero, and unused component entries must retain defaults. Component energy
uses the integrated heat-capacity polynomial relative to its reference temperature.
The surface cap limits the pure-component vapor-state temperature before Raoult
partitioning; it does not implement physical mixture boiling. Models, coefficients
and bounds enter restart identity. Default constant-property contexts are preserved.
See [component property evidence](validation/spray-mixture-properties-20260916.md).

Initial and scheduled parcels use `mixture_fractions` by default. To override a
port, set `port_mixture_fractions(component,port)`, a `(16,16)` array. For example,
`port_mixture_fractions(1:2,1)=.2,.8` assigns the first port a 20/80 mass blend.
Each nonzero column must be nonnegative and sum to one; an all-zero column inherits
the shared composition. Rows beyond `mixture_count` and columns beyond
`injection_port_count` must stay zero. Components may be absent from individual
ports; every port uses the same ordered component property table. Keep a valid
shared composition even when overriding all active ports.

Both initial and scheduled injections use each port's selected composition for
density, energy, component inventories and rate-derived velocity. Restart identity
includes active port compositions. See [port-composition evidence](validation/spray-port-mixture-20260916.md).
Existing size laws,
nozzle geometry, per-port temperatures and shared/independent schedules remain
available. Injection density and energy use the actual mixture at each port's
injection temperature; rate-derived velocity uses that density. New component
masses append alongside each injected parcel, with species-resolved mass ledgers.
Per-port pure-liquid property overrides cannot be combined with mixture mode. Prescribed mixture viscosity/surface tension permit representative
ETAB/KHRT, KH shed children and computational splitting.

Parcel CSVs add `mass_<species>` columns in kg per physical droplet. Their sum is
`mass`; `diameter` and Sauter mean diameter use the evolving mixture density.
The ordered component table and initial fractions enter restart compatibility.
Changed composition/properties reject restart, and checkpoints retain every
parcel's evolved component masses. The terminal/trace update changes the
mixture evolution tag to `mixed-liquid-components-v2-terminal`; earlier mixture
checkpoints reject rather than silently changing their continuation policy.
See [terminal-state evidence](validation/spray-mixture-terminal-20260916.md).

[mixture.nml](../cases/spray_3d/mixture.nml) is a small synthetic example. The
illustrative H2O/H2O2 properties demonstrate software behavior, not a validated
liquid blend. See [CLI/injection evidence](validation/spray-mixture-cli-20260916.md).

### Mixed-liquid wall-impact event API

`resolve_mixture_wall_impact` partitions an incident mixture among the retained
parent, splash children and a deposited film patch using the existing six-regime
impact kernel. Parent and child component masses are kg per physical drop;
`film_mass` is kg per complete film patch. All outgoing material inherits the
incident composition and temperature. The kernel evaluates mixture density and
energy at that temperature and returns wall momentum/energy exchange. Failed
classification or child/film construction preserves all supplied state arrays.

The caller supplies the dimensionless impact parameters, including wall-to-boiling
temperature ratio. The planar dispatcher below evaluates that ratio. Grid/CLI coupling remains to be connected. The checkpoint API below retains
mixture-film composition. The CLI supports this mixed-liquid wall path. There is no entrainment
from an existing film. See [event evidence](validation/spray-mixture-wall-event-20260916.md).

### Mixed-liquid film exchange API

`advance_mixture_wall_film_exchange` accepts component evaporation amounts in kg,
component vapor enthalpies in J/kg and gas/wall heat amounts in J into a stationary
film patch. It inverts the remaining mixture energy for temperature and updates
height at fixed footprint using the new density. `component_mass` stores total kg
for the patch. Returned gas changes include individual vapor species; `wall_energy`
is the opposite of supplied wall heat. Momentum exchange is zero for this stationary
film operation.

Complete evaporation and an evaporation-driven height below `minimum_height` return
all remaining species and energy to gas. The full inventory convention also applies
when the final temperature cannot be defined because no liquid remains. Invalid
component loss or temperature inversion leaves the patch and composition untouched.
This supplied-exchange API is used by the frozen-rate update below. Mixed films
are not yet connected to the grid driver. See
[film exchange evidence](validation/spray-mixture-film-exchange-20260916.md).

### Mixed-liquid film rates and frozen update

`mixture_wall_film_rates` evaluates Raoult surface fractions and diffusion-weighted
component evaporation through the existing stationary-film correlation. Gas-side
properties use the wall half-composition rule and
`Tskin = Tfilm + 0.5*max(0,Tgas-Tfilm)`. Gas and wall heat rates use the film thickness
and supplied liquid conductivity. The temperature-rate estimate subtracts each
component's evaporation rate times its own vapor-minus-liquid enthalpy, divided by
the film mass and mixture heat capacity.

`advance_mixture_wall_film_frozen` holds these rates fixed over a supplied interval
and calls the conservative exchange API. Actively evaporating components whose
current or predicted residual mass is at most `1e-12` of film mass are returned
entirely to gas. Ambient-suppressed components stay in the film; a zero interval
changes nothing. The grid depletion limiter excludes already trace components.
Bulk component overdraw or an
inadmissible remaining temperature rejects; a grid driver must reduce the step in those
cases. It retains the full-inventory minimum-height removal policy. Composition
transport, thermal updates and rates are now available as APIs; grid/CLI feedback
and application restart wiring remain open. See
[rate evidence](validation/spray-mixture-film-rates-20260916.md).

### Mixed-liquid planar wall dispatch

`mixture_impact_at_planar_wall` finds the earliest crossed box face, calculates
impact parameters with the current mixture density and prescribed mixture viscosity
and surface tension, and dispatches the existing impact event. Outgoing parcel
component arrays use kg per drop; deposition components use whole-patch kg.
No crossing leaves the supplied state unchanged; failed dispatch rolls back.

`mixture_impact_boiling_temperature` uses the mass-weighted pure-component boiling
temperatures for the wall-temperature ratio, following pinned `boilT_mix` and
AhamedSplash. Component temperatures use the existing Clausius or bounded Antoine
inversion at local pressure. This empirical impact temperature is not a Raoult
bubble-point calculation or a physical boiling-rate model. Existing film geometry
supplies local thickness; film entrainment remains absent. See
[crossing evidence](validation/spray-mixture-wall-crossing-20260916.md).

### Mixed-film checkpoint API

The checkpoint reader/writer accept `film_component_mass(component,patch)` in kg
per complete film patch. It accompanies `component_mass`, `film_state=.true.` and
`wall_films`; the format adds `_FMIX_V1` and a component row after each film record.
The number and order of components match the particle component matrix. Callers
must bind the ordered component property table and film settings in the context.

A mixed checkpoint with film state requires this matrix, including a zero-column
matrix when there are no films. Shape mismatch, negative/nonfinite component mass
or failure to sum to film mass rejects. Reading validates the complete candidate
before changing gas, particles, films or either component matrix. Existing pure-film
and particle-only mixed formats retain their layout. The CLI retains this mixed-film checkpoint state.
See [checkpoint evidence](validation/spray-mixture-film-checkpoint-20260916.md).

### Mixed-film grid feedback API

`advance_mixture_film_coupling_3d` advances stationary film patches and their
component matrix against a regular gas grid. It samples gas properties at each
patch, limits substeps by component depletion and temperature change, and deposits
component mass and energy with the particle interpolation weights. Failed gas
recovery halves the attempted interval; failure of the complete update preserves
gas, temperatures, films and component masses. Optional `start_times` delays each
patch until its birth time within the interval.

The returned boundary ledger records heat entering the simulated gas/film system
from walls. Fully evaporated patches and corresponding matrix columns are removed
together. `mixture_system_totals` accepts optional `wall_films` and `film_mass` to
include stationary-film mass, species and internal energy. This API is connected to particle impact and reactive/CLI stepping. See [grid evidence](validation/spray-mixture-film-grid-20260916.md).

### Combined mixed-particle and wall-film update

`advance_mixture_coupling_3d` accepts optional `wall_options`, `wall_films` and
`film_component_mass`. The two film arrays must accompany each other. Enabled wall
impact requires these arrays, even when initially empty. It processes the accepted
particle flight's first wall crossing before aerodynamic breakup; a hit skips
breakup for that substep. Splash children and deposited patches inherit the incident
composition and start at the accepted substep's end.

After particle flight, stationary films exchange mass and energy with the grid for
their remaining interval. Particle and film substeps share the configured budget.
Capacity failure, inadmissible impact or rejected film feedback rolls back the entire
particle/gas/film update. Boundary ledgers combine wall impulse/energy and wall heat.
Reactive-driver and CLI integration are available below. See
[combined evidence](validation/spray-mixture-wall-coupled-20260916.md).

### Reactive mixed-wall integration

`advance_spray_reactive_3d` accepts `film_component_mass` alongside mixed-particle
composition and `wall_films`. Both spray half-steps now carry wall impact and mixed
film evolution through the same candidate state as hydro, molecular transport,
LES and chemistry. The complete state commits only after every stage succeeds;
a late chemistry failure restores both component matrices with gas, particles and
films. Pure-material property tables still cannot accompany mixture tables.

The existing prohibition on moving particle-impact walls remains. CLI input/output
and checkpoint wiring retain both particle and film compositions. See
[reactive evidence](validation/spray-mixture-wall-reactive-20260916.md).

Mixed-film evolution now uses `mixture_film_evolution_id='mixed-film-v2-trace'`
in CLI checkpoint context. Earlier mixed-wall checkpoints reject because their
trace-depletion evolution differs. The storage layout, pure-wall and particle-only
mixture contexts are unchanged. API callers must include the evolution ID in their
own context. See [trace evidence](validation/spray-mixture-film-trace-20260916.md).

### Mixed-wall CLI and outputs

Combine mixture inputs with the existing `wall_impact_faces`, `surface_tension`,
`liquid_viscosity`, `wall_film_evolution` and `wall_liquid_conductivity` settings.
Stationary film heat conductivity is prescribed; this is not a calibrated mixture
conductivity model. The `.films.csv` output adds `mass_<species>` columns in kg per
whole patch, while parcel component columns remain kg per physical drop. History
and conservation ledgers include gas, particles and mixed films.

Restart retains both component matrices and validates film volume using the current
film composition and temperature. [mixture_wall.nml](../cases/spray_3d/mixture_wall.nml)
is a small synthetic deposition/evaporation example. See
[CLI evidence](validation/spray-mixture-wall-cli-20260916.md).

## KH-RT representative-drop breakup

Set `breakup_model='khrt'` with positive `surface_tension` and
`liquid_viscosity`. `khrt.nml` uses the pinned PelePhysics defaults
`khrt_b0=.61`, `khrt_b1=7`, `khrt_c3=1`. KH waves relax the representative
radius; RT instability uses a persistent waiting time, initialized to -1
and enabled by the first KH event. The Reynolds number uses film density and
viscosity; the wave equations use bulk gas density and relative velocity.
Mass reduction increases multiplicity at fixed represented mass, velocity and
temperature. The minimum diameter and transactional failure rules still apply.

The default implements the upstream `make_new_drops=false` branch. Optional
KH shed-mass children are described below. Computational splitting partitions
identical representative droplets. Surface energy is omitted.
Zero RT acceleration resets an active timer to avoid singular expressions.
Exact upstream timestep parity and physical validation have not been established.

Parcel CSVs include `rt_time`. KH-RT checkpoints use V7/V8 (periodic/physical)
or V9/V10 with computational splitting, retaining the RT timer and binding all
three coefficients into the restart context. Existing V1--V6 formats are unchanged.
See [source provenance](../references/spray_khrt.json) and
[focused evidence](validation/khrt-20260915.md).


### KH shed-mass children

Set `khrt_children=.true.` and positive `khrt_reference_diameter_cubed`
with `breakup_model='khrt'`. The reference quantity has units of m³ and is
`mean_over_jets(N_injection * mean_diameter**3)`, matching the upstream
`avg_inject_mass` normalization despite that variable's name. Supply it explicitly;
it is neither liquid mass in kg nor an automatically inferred injection dose.

The model accumulates `N_initial * (r_old**3-r_new**3)` during KH events.
This inventory remains included in the parent until release. The reference's
3% jet threshold and 90% parent threshold select partial or full release;
partial release also requires a resolved child radius and sufficient child count.
The volume denominator's numerical epsilon is converted from cm³ to m³.
The configured `minimum_diameter` supplies the radius floor, rather than the
reference's density-dependent minimum-mass rule.

Each release creates four equal-weight children in opposite transverse pairs,
with speed `0.188 * omega_KH * lambda_KH`. They inherit the parent's position,
temperature and liquid material, retain `parent_id`, and start fresh breakup states.
The transverse plane uses the parent velocity (relative gas velocity if stationary).
Gas energy decreases by the added transverse kinetic energy, closing total energy;
the paired directions preserve mean momentum. This deterministic quadrature is a
PeleF choice and does not reproduce the upstream random azimuth/count algorithm.
Creation occurs at the accepted parcel substep's end; children advance through
any remaining interval. Failed feedback, capacity checks or later chemistry roll
back the entire coupled update.

CSV fields `khrt_initial_count` and `khrt_shed_volume` expose the accumulator;
checkpoint suffix `_SHED_V1` retains it alongside RT time and parent lineage.
The enable flag and reference quantity participate in restart compatibility.
Stationary wall impact and computational splitting are supported. For splitting,
both KH accumulators are partitioned with the represented drop count, so a split
does not duplicate pending shed mass. The final piece carries rounding remainders.
The jet release threshold still applies per parcel against the prescribed global
reference; changing computational multiplicity can therefore change release timing.
When a substep crosses an enabled wall, its incident parcel undergoes wall impact
before any KH release. A hit skips aerodynamic breakup for that substep; otherwise
KH children are born at its end and can impact during their remaining flight.
This prevents the same incident mass from being released and deposited twice.
Splash children start with empty KH inventories; deposition consumes the parent's
inventory, thermal breakup resets it, and a rebound retains it. Existing wall
correlation rejection rules still apply, including unsupported wet-splash states.
See [combined wall/KH evidence](validation/khrt-wall-20260916.md).
RT still changes the representative size; it does not itself create children.
Collision, general turbulent dispersion, primary atomization and physical
qualification remain open. See [focused evidence](validation/khrt-children-20260916.md).

### Computational parcel splitting

For ETAB or KH-RT runs, set `maximum_multiplicity` to a positive upper bound on represented
droplets per computational parcel. Zero (the default) disables this partition.
`maximum_parcels` bounds allocated particle storage, up to 1,000,000.
After a coupling interval, overweight live parcels are partitioned into equal
weights, with the final child carrying the rounding remainder. Droplet mass,
position, velocity, temperature, deformation and RT timer are copied. KH initial
count and pending shed volume are divided among the pieces; aggregate mass,
momentum and energy are conserved. This increases computational parcel count
without introducing a transverse velocity or dispersion model.

Existing IDs remain fixed. New IDs append in parent order; each new child stores
the ID of the parcel it was split from. Injection roots have `parent_id=0` and
later injections continue after the last allocated ID. Restart validates parent
links and counts roots against injection events, independently of child count.
The complete split plan is checked against storage and ID limits before arrays
are replaced; failure rolls back the coupled gas/particle update.

ETAB splitting-enabled checkpoints use V5/V6, storing exact integer parent IDs in
addition to the ETAB and optional boundary state. The multiplicity/storage
limits are included in the context; parcel CSVs gain `parent_id`. For example,
add `maximum_multiplicity=1e7, maximum_parcels=1000` to `etab.nml`.
See [splitting/restart evidence](validation/spray-splitting-20260915.md) and
[KH shedding partition evidence](validation/khrt-splitting-20260916.md).


## Gas reservoir inlet

Use `boundary_faces='inflow','outflow','periodic','periodic','periodic','periodic'`
with `inflow_temperature(1)`, `inflow_pressure(1)` and `inflow_velocity(:,1)`
for an x-low gas reservoir. The six-face ordering is unchanged. The defaults are
1200 K, 101325 Pa and zero velocity; active inlet values must be finite and
thermodynamically valid. The state is prescribed in ghost cells for hydro,
molecular transport and SGS. This is not a characteristic nonreflecting inlet
or a specified numerical mass-flux boundary.

Inlet composition defaults to the initial gas mixture. Override it with
`inflow_composition_count`, `inflow_composition_species(32)` and
`inflow_composition_mole_fractions(32)`, using named species of the loaded
mechanism. A shared inlet mixture applies to all inlet faces; thermodynamic
state and velocity may differ by face. `inflow.nml` exercises an H2/O2/N2
reservoir feeding initially pure N2 with evaporating water parcels, molecular
transport and Cs=0.17.

The timestep is capped by the inlet wave speed and, when enabled, its molecular
diffusion limit. Gas convective/diffusive/SGS exchange joins the surface ledger;
all active inlet settings and resolved composition bind checkpoint context.
Particles crossing out through an inlet leave the system and are accounted as
particle outflow. A gas inlet does not create parcels; parcel injection remains
controlled by the existing separate source schedule.

See [inlet restart and conservation evidence](validation/spray-inflow-20260915.md).


## Isothermal and tangentially moving walls

`wall_temperature(6)` specifies kelvin for `isothermal_wall` faces.
`wall_velocity(3,6)` specifies a tangential velocity on `no_slip_wall` or
`isothermal_wall` faces; its normal component must be zero. Defaults are
300 K and zero velocity. Heated or moving walls require molecular transport.
`thermal_moving_wall.nml` supplies a small executable case.

Molecular wall heat transfer and viscous work use the existing regular-grid
transport operator. SGS uses the same reflected wall velocity and applies a
wall-to-cell thermal gradient over half a cell. Thermal gradients enter the
SGS timestep limit; the molecular wall limit also caps the application step.
All wall fluxes enter the accumulated momentum/energy ledger. Active wall
parameters bind restart context, so changed wall motion or temperature rejects
an incompatible restart.

Parcel impact remains specular: normal velocity reverses, tangential velocity
and temperature remain unchanged. Wall heat exchange with droplets, splash,
adhesion, roughness and wall-function calibration remain open. The gas wall
motion is tangential motion at a fixed boundary, not moving geometry.
See [focused checks](validation/spray-thermal-wall-20260915.md).


## Uniform external acceleration

`acceleration=0,0,-9.81` applies a prescribed constant acceleration in m/s² to
both gas and droplets; the default is zero. See `acceleration.nml`. Two exact
half-step velocity kicks surround the coupled operators. Gas total energy
receives the corresponding kinetic work; droplet temperature is unchanged.
The midpoint particle drift is exact for constant acceleration without drag.
The overall evaporating, coupled method retains its existing accuracy limits.

External impulse and work are separate from gas/liquid feedback and physical
boundary exchange. History gains `source_px`, `source_py`, `source_pz` and
`source_energy`; conservation residuals subtract these accumulated sources.
Checkpoints with acceleration append `_FORCE_V1` to their existing format tag
and store a source ledger. The acceleration vector binds restart context.
Zero-acceleration files retain their previous format. Step rejection rolls back
both phases and external-source increments. The timestep includes an
acceleration displacement limit, and kicks check thermal-state recovery.

This is a uniform body force, not a hydrostatic well-balanced discretization,
spatial gravity field, buoyancy closure or particle-only external force model.
Potential energy is not included in reported energy; external work accounts for
its effect. See [focused evidence](validation/spray-acceleration-20260915.md).


## Larger selected mechanisms

The importer, generator, selected input arrays, reactive state/transport checks
and MPI species-comparison storage now share a 128-species capacity. Runtime
Fortran limits use `constants_mod::mechanism_species_capacity`. Reaction and
thermo families remain restricted to the documented supported forms.
The capacity is a storage/input bound, not a performance qualification.

A 53-species, 325-reaction GRI-Mech case now runs through the existing selected
CVODE 0D and regular 3D spray executables. The 0D trajectory agrees with an
independent Cantera reactor, and the coupled chemistry/transport/spray case
restarts exactly. See [commands and evidence](validation/species-capacity-20260915.md).
MPI and AMR with more than 32 species have not been newly qualified.


## SRI pressure-dependent reactions

Selected mechanism ingestion accepts Cantera `SRI` falloff data with required
`A`, `B`, `C` and optional `D=1`, `E=0`. The normalized bundle uses
`type="sri-falloff"` and an `sri` coefficient mapping. SRI is distinct from
Troe and cannot specify both models. A and C must be nonnegative, D positive,
and A and C cannot both vanish. All coefficients must be finite.

Fortran evaluates the five-parameter broadening factor using log-sum-exp and
provides the analytical concentration derivative used by the existing kinetic
Jacobian. AMR/EB checkpoint mechanism records include the five SRI coefficients
for SRI reactions; MPI mechanism comparisons include them too. Existing
non-SRI checkpoint records keep their layout. This extends reaction kinetics,
not the liquid property model. See [focused evidence](validation/sri-20260915.md).


## Chemically activated reactions

Cantera `type: chemically-activated` reactions now use the same supported
Lindemann, Troe, SRI or Tsang broadening data as falloff reactions. Their rate is
`k0 F / (1 + Pr)`, with `Pr = k0 [M] / kinfinity`; the runtime supplies the
corresponding analytic concentration derivative. This is distinct from the
ordinary falloff prefactor. Normalized bundles use `chemically-activated`
or `sri-chemically-activated` as their reaction type.

Reaction kinds identify the different pressure laws in MPI mechanism comparison
and checkpoint records. SRI coefficients persist for both pressure laws, and a
changed coefficient rejects AMR restart. The existing reduced-pressure floor
also applies to the broadening calculation; below it the broadening derivative
is zero for the activated branch. See [focused evidence](validation/activated-rates-20260915.md).


## Tsang broadening

The selected importer accepts Cantera `Tsang: {A: ..., B: ...}` data for both
falloff and chemically activated reactions. B defaults to zero. Normalized
bundles use `tsang-falloff` or `tsang-chemically-activated` and a `tsang`
coefficient mapping. The center factor is `A+B*T`; the rest of the broadening
and its concentration derivative use the existing Troe expressions. Both
coefficients must be finite and the center must be positive at the evaluated
temperature. A nonpositive center rejects the rate calculation.

Tsang has separate reaction kinds and coefficient storage. AMR/EB checkpoint
records and MPI mechanism comparisons include both coefficients. Legacy
non-Tsang checkpoint layouts are preserved. See
[focused evidence](validation/tsang-20260915.md).


## PLOG pressure interpolation

The importer accepts Cantera `pressure-dependent-Arrhenius` reactions and emits
normalized `type="plog"` with ordered `plog` entries containing P (Pa), A, b and
Ea (J/kmol). Up to 64 entries per reaction are stored. Equal-pressure entries
remain separate so their temperature-dependent rates can be summed at runtime.
Signed A coefficients are supported when each pressure group's total is positive.
Nonpositive group totals or unrepresentable rates reject the calculation.

The runtime sums each group using scaled signed exponentials, then interpolates
log(rate) against log(pressure). Outside the supplied pressure interval it uses
the nearest endpoint rate. Pressure is R*T*sum(concentrations); its analytic
concentration derivative feeds the existing Jacobian. Internal knot derivatives
use the right interval, and endpoint branches have zero pressure derivative;
these derivatives need not be continuous at knots.

PLOG count and coefficients participate in MPI mechanism comparison and
AMR/EB checkpoint records. Existing non-PLOG checkpoint records retain their
layout. See [focused evidence](validation/plog-20260915.md).

## Chebyshev temperature/pressure rates

Selected bundles accept Chebyshev coefficient matrices up to 8 by 8, with
positive, ordered temperature and pressure fit bounds. Coefficients come from
Cantera's normalized rate API in m/kmol/s units. The runtime evaluates the
polynomial in inverse temperature and log pressure, with pressure R*T*sum(C),
and supplies its analytic derivative with respect to every concentration.
It evaluates the polynomial outside the fit interval when representable;
this extrapolation does not extend the physical validity of a mechanism fit.
Zero pressure, malformed tables and unrepresentable logarithmic rates fail.

The matrix, shape and bounds participate in MPI mechanism comparisons and
AMR/EB checkpoint contexts. Existing reaction kinds retain their checkpoint
layouts. See [the focused Chebyshev evidence](validation/chebyshev-20260915.md).

## Lu n-dodecane gas mechanism

`prepare_spray_reference.py --fuel dodecane` prepares the pinned 53-species,
268-reaction Lu mechanism using the same source/license checks as the default
methanol reference. Configure the selected bundle with its generated JSON and
original YAML. The importer and runtime preserve finite Troe A values outside
(0,1), including A=1, and signed T2; T1/T3 remain positive and the evaluated
center must be positive and representable.

The selected CVODE path now evaluates negative Newton trial concentrations
through a nonnegative extension and matching Jacobian columns. Accepted states
retain the previous composition, conservation and CVODE constraints. This
allows trace species to approach zero during stiff ignition without rejecting
the internal Newton iteration prematurely.

A phi=0.5, 1100 K, 4 MPa constant-volume ignition case reaches 2478 K at 0.6 ms
and agrees with Cantera. A separate short coupled case exercises this gas
mechanism with water parcels and exact restart; it does not supply a liquid
n-dodecane property model. See [commands and evidence](validation/dodecane-20260915.md).

## Antoine saturation-pressure fits

Set `vapor_pressure_model='antoine'`, `antoine_coefficients=A,B,C,D`, and
`antoine_temperature_range=Tmin,Tmax` to use
`psat = D * 10**(A - B/(T+C))`. T, B and C use kelvin; D scales the resulting
pressure to Pa. For coefficients tabulated in bar, D=100000. The convention
matches the pinned PelePhysics `SprayPropertiesPeleMP.H` implementation.

The default `vapor_pressure_model='clausius'` retains the existing relation
and requires the Antoine arrays to be zero. Antoine inputs require finite
coefficients, positive B/D, ordered positive bounds and T+C>0 throughout the
interval. Evaluation outside the fit range or at psat>=ambient pressure fails
transactionally. This does not implement boiling, supercritical evaporation,
multicomponent liquid properties. Temperature-dependent heat capacity is configured separately below.

`cases/spray_3d/water_antoine.nml` supplies the NIST water fit for 273--303 K
and starts droplets at 290 K. The coefficient source and units are recorded in
[`references/spray_antoine.json`](../references/spray_antoine.json); density/cp
remain illustrative constants. Changing either coefficients or bounds rejects
a checkpoint restart. See [focused evidence](validation/antoine-20260915.md).

## Temperature-dependent liquid density

Use `liquid_density_model='polynomial'`, `liquid_density_coefficients=a,b,c,d`
and `liquid_density_temperature_range=Tmin,Tmax` for
`rho(T)=a+T*(b+T*(c+T*d))` in kg/m3, with T in kelvin. This is the pinned
PelePhysics `MPLiqProps::rho_i` convention. The coefficients must already use SI
units. `liquid_density` remains the constant-model input and initialization
reference; the polynomial supplies the actual density when selected.

Injection converts specified diameter to mass at the injection temperature.
Drag, frozen-film evaporation and ETAB/KH-RT use the current droplet density.
Evaporation freezes density over each substep, while the resulting diameter
uses the updated temperature. Heating can therefore change volume without
creating liquid mass. The formation-energy reference remains unchanged; liquid p/rho and surface-energy corrections are
still neglected in this dilute point-particle model.

Coefficients and bounds participate in restart identity. Evaluating outside
the interval or obtaining density outside 1e-3--1e5 kg/m3 rejects the update;
a polynomial is not guaranteed positive over its entire interval by setup.
The default constant-density path retains its checkpoint context.
`cases/spray_3d/water_density.nml` uses an illustrative linear polynomial to
exercise the implementation, not calibrated water properties. See
[focused evidence](validation/liquid-density-20260915.md).

## Temperature-dependent liquid viscosity for breakup

`liquid_viscosity_model='inverse-temperature'` selects
`mu(T)=a+b/T+c/T**2+d/T**3` with `liquid_viscosity_coefficients=a,b,c,d` and
`liquid_viscosity_temperature_range=Tmin,Tmax`. Use kelvin and coefficients
that return Pa s. The convention matches pinned PelePhysics `MPLiqProps::mu_i`.
ETAB damping and KH-RT's Ohnesorge number use this viscosity at the current
droplet temperature. The model requires an enabled breakup model. It can be
combined with temperature-dependent density.

The default remains `liquid_viscosity_model='constant'`, using
`liquid_viscosity`; its new coefficient/range arrays must be zero. With the
inverse-temperature model, the scalar is retained in input identity but the
fit determines viscosity. Evaluations outside the temperature range or outside
1e-10--1e3 Pa s fail transactionally. Fit coefficients and bounds participate
in restart identity. No positivity proof across the whole fit range is implied.

`cases/spray_3d/water_viscosity.nml` exercises combined density/viscosity and
KH-RT using illustrative coefficients, not a calibrated water model. See
[focused evidence](validation/liquid-viscosity-20260915.md).

## Near-boiling vapor-state temperature cap

`boiling_model='surface-cap'` with `critical_temperature=Tc` enables the pinned
PelePhysics vapor-state convention `Ts=min(Tdrop,0.99*min(Tboil(P),Tc))`.
The pressure-dependent boiling temperature comes from the existing constant
reference latent heat and Clausius--Clapeyron equation. The cap affects the
saturation-pressure evaluation; the droplet temperature, energy accounting,
film transport temperature, and density/viscosity fit bounds retain their
existing meanings.

This option currently requires `vapor_pressure_model='clausius'` and Tc>Tref.
It retains the default `boiling_model='reject'` behavior when absent. A
nonpositive inverse boiling temperature is rejected. The critical temperature
and selected model are part of restart identity. The example
`cases/spray_3d/water_surface_cap.nml` starts water droplets at 400 K.

This is a near-boiling numerical vapor-state treatment, not a heat-limited
boiling-rate, supercritical-fluid or bubble-growth model. The upstream Watson
latent-property correlation is not added here; the supplied reference latent
heat remains authoritative. See [source and evidence](validation/surface-cap-20260915.md).

## Multiple injection ports

Set `injection_port_count=1..16`, then provide `injection_port_centers(3,n)`,
`injection_port_axes(3,n)` and `injection_port_fractions(n)`. In a namelist,
list x/y/z for the first port, then x/y/z for the next. Fractions must be
positive and sum to one (roundoff is normalized). Unused entries must be zero.
Each active center must lie inside the domain and each axis must be nonzero.

`initial_liquid_mass` and `mass_flow_rate` are totals over all ports. Each port
receives its fraction of the event mass and creates `parcel_count` parcels.
IDs advance in port order, preserving deterministic restart and parent counts.
By default the schedule, liquid, temperature, speed, cone angle and diameter
distribution are shared. Independent schedules are available below; distinct
fuels per port and resolved nozzle flow remain outside this path. `maximum_parcels` bounds the combined allocation.

The default count zero uses the original single `injection_center`/`injection_axis`
path and checkpoint context. With active ports, their centers/axes replace the
single geometry; all port data participate in restart identity. See
`cases/spray_3d/water_ports.nml` and [focused evidence](validation/injection-ports-20260915.md).

## Independent port schedules

Set `independent_port_schedules=.true.` with active injection ports, then supply
`port_mass_flow_rates`, `port_injection_starts`, `port_injection_ends` and
`port_injection_intervals` in port order. The global `mass_flow_rate` must be
zero. Each active port needs a positive interval and a nonempty time window.
A zero rate disables scheduled injection at that port. Unused entries are zero.
The aggregate input is limited to one million scheduled events.

Each port injects `rate*min(interval,end-event_time)` and creates `parcel_count`
parcels per event. Events use start+n*interval and occur strictly before end.
The solver processes all due events before stepping and breaks exact timestamp
ties by port number. Initial liquid mass still follows `injection_port_fractions`;
those fractions do not multiply the independent scheduled rates. Material remains shared. Speed, temperature, diameter and cone angle can be
overridden per port as described below.

Restart preserves the aggregate event count and reconstructs the per-port
counts by replaying the deterministic schedule. It checks injected mass,
parent counts, future events and missed events before continuing. Schedule
arrays are included in checkpoint identity. Replay is bounded by the input's
one-million-event limit; large restart performance is not qualified.
See `cases/spray_3d/water_independent_ports.nml` and
[focused evidence](validation/independent-port-schedules-20260915.md).

## Per-port injection conditions

`injection_port_speeds`, `injection_port_diameters`,
`injection_port_temperatures` and `injection_port_cone_angles` override the
corresponding global values for active ports. An entry of -1 (the default)
inherits its global value; unused entries must remain -1. Units are m/s, m, K
and radians. The existing range/material checks apply to each resolved port.

Initial injection and both shared and independent schedules use the resolved
conditions. Different temperatures and diameters produce the corresponding
single-drop mass and multiplicity, while each event retains its specified
mass allocation. Momentum and energy ledgers use each generated parcel's own
velocity and temperature. By default the liquid, distribution shape and truncation bounds
remain shared. For a distributed diameter input, the per-port diameter replaces
the same scale parameter as `drop_diameter`.

Resolved override values participate in restart identity. Without overrides,
the preceding port checkpoint context remains unchanged. See
`cases/spray_3d/water_port_states.nml` and
[focused evidence](validation/injection-port-states-20260915.md).


## Different pure liquids by injection port

`port_vapor_names(16)`, `port_liquid_densities(16)`, `port_liquid_cps(16)`,
`port_reference_temperatures(16)`, `port_reference_pressures(16)` and
`port_latent_heats(16)` override each active port's liquid. Units are SI;
blank names and numeric -1 inherit the global values. Unused slots must keep
those defaults. Every vapor name must exist in the configured gas mechanism;
formation energy is initialized separately from that species' thermodynamics.
Specifying any override enables material IDs equal to the port number.

Initial injection and both schedule modes retain this identity. Droplet mass,
heat transfer, evaporation, gas-species feedback, boundary budgets and system
energy use the selected liquid. Computational splitting copies the identity.
Parcel CSV adds `material_id`; opt-in checkpoint records retain it and reject
changed liquid inputs. Default single-liquid checkpoint/output formats remain.

Each port can also select its density and saturation-pressure models and
near-boiling cap, as described below. Breakup viscosity and surface tension can also be overridden per port;
model selection and the remaining breakup controls stay shared. The low-level material
pool API accepts individual `spray_liquid` records; totals require validated
parcel IDs and materials, as supplied by injection/coupling or checked after
application restart. This is coexisting pure-liquid droplets, not a liquid
mixture with composition-dependent evaporation.

`cases/spray_3d/port_materials.nml` uses H2O and H2O2 in inert nitrogen with
chemistry disabled. Its second-liquid properties are illustrative integration
inputs, not a calibrated peroxide model.
See [focused evidence](validation/port-materials-20260915.md).


### Property models per port

The following arrays select individual liquid models (up to 16 ports):

| Input | Shape | Default |
| --- | --- | --- |
| `port_density_models` | 16 | `constant` |
| `port_density_coefficients` | 4,16 | 0 |
| `port_density_temperature_ranges` | 2,16 | 0 |
| `port_vapor_pressure_models` | 16 | `clausius` |
| `port_antoine_coefficients` | 4,16 | 0 |
| `port_antoine_temperature_ranges` | 2,16 | 0 |
| `port_boiling_models` | 16 | `reject` |
| `port_critical_temperatures` | 16 | 0 |

Select `polynomial`, `antoine` or `surface-cap` as for the corresponding global
model. Coefficients, SI units, bounds and rejection rules are identical.
For example, `port_density_coefficients(:,2)=1150,-.5,0,0` and
`port_density_temperature_ranges(:,2)=280,500` supply the second port's
illustrative linear density fit after selecting its `polynomial` model.

Any nondefault port model input activates per-port material identity even if
all vapor names inherit the global name. Port model fields have their own
explicit defaults; they do not inherit global fitted coefficients or model
selectors. This prevents a fit for one global liquid from being applied to a
different port liquid implicitly. Unused slots retain the defaults. A cap still
requires Clausius pressure and a critical temperature above the liquid reference
temperature; Antoine and cap cannot be combined for the same material.
All port model fields participate in restart identity when enabled. Default
material contexts remain compatible. See
[focused evidence](validation/port-property-models-20260915.md).


### Breakup properties per port

`port_surface_tensions(16)` and `port_liquid_viscosities(16)` override surface
tension (N/m) and dynamic viscosity (Pa s). Numeric -1 inherits the global
value. `port_viscosity_models(16)` defaults to blank, inheriting the complete
global viscosity model and fit. An explicit `constant` or `inverse-temperature`
selector instead uses that port's `port_viscosity_coefficients(4,16)` and
`port_viscosity_temperature_ranges(2,16)`, both defaulting to zero. Coefficients
without an explicit selector are rejected. Inverse-temperature units, validity
bounds and the evaluation formula match the global viscosity model.

These inputs require active ports and ETAB or KH-RT breakup. Keep valid global
breakup settings as inheritance defaults. Unused port fields must retain their
defaults; injection temperature must lie in the chosen viscosity fit's domain.
Later out-of-domain candidates retain the existing whole-step rollback behavior.
Any override enables material identity, including when vapor species are shared.
Both coupling half steps select breakup properties by that identity. Computational
splitting preserves identity; model selection, calibration constants and global
storage/substep controls remain shared in application inputs. Resolved properties
are included in restart identity. With no overrides, preceding formats remain.

The core coupling API accepts an optional `breakup_options` array alongside the
material pool, with one entry per material and the same breakup model as global
options. Only its breakup fields are used; global options still control transfer,
substeps and computational splitting. See
[focused evidence](validation/port-breakup-properties-20260916.md).


## Temperature-dependent liquid heat capacity

Use `liquid_cp_model='polynomial'`, `liquid_cp_coefficients=a,b,c,d` and
`liquid_cp_temperature_range=Tmin,Tmax` for
`cp(T)=a+b*T+c*T**2+d*T**3` in J/(kg K), with kelvin temperatures. The default
`constant` uses `liquid_cp` and requires zero coefficient/range arrays.
The polynomial must remain within [1,1e8] J/(kg K) over its complete interval;
validation checks endpoints and all interior extrema. Tref must lie inside
that interval. Coefficients are bounded by 1e10 in absolute value; temperatures
are within [1,1e5] K. Injection and subsequent live droplets cannot leave the
fit interval. Out-of-range update candidates retain whole-step rollback.

Liquid energy uses the exact polynomial integral from Tref, anchored to the
existing gas enthalpy minus reference latent heat. Both the heating/evaporation
update and its substep estimate use cp at the substep's initial droplet
temperature. This is a frozen-property time approximation; it does not exactly
integrate the nonlinear temperature ODE. Discrete gas feedback uses the change
in integrated liquid energy, preserving all-phase total energy. Clausius
saturation pressure still uses the prescribed reference latent heat; this
change does not add a variable-latent-heat saturation model.

For separate port materials use `port_cp_models(16)`,
`port_cp_coefficients(4,16)` and `port_cp_temperature_ranges(2,16)`. These
independently default to constant/zero, like other per-port property models;
constant cp inherits the existing per-port/global scalar cp input. Any nondefault
port cp field enables material identity. Unused slots retain defaults. All fits
and intervals enter checkpoint identity; default constant-cp formats remain.

The tests use an illustrative `cp=4000+4*T` fit over 300--600 K, with no claim
of measured liquid properties. See [focused evidence](validation/liquid-heat-capacity-20260916.md).


## Ahamed wall-impact correlations

`spray_wall_impact_mod` ports three routines from the pinned PelePhysics
`AhamedSplash.H`: impact-regime selection, five splash-velocity components and
four azimuthal mass fractions. Regimes cover rebound, deposition, dry splash,
thermal breakup, spreading and wet splash. The API accepts the dimensionless
impact parameter `Kv`, wall/boiling temperature ratio, impact angle (30--90
degrees in radians), critical impact parameter and film-thickness/drop-diameter
ratio. Velocity inputs/outputs use m/s. Azimuthal fractions are fractions of
incident mass, summing to the supplied splashed fraction when accepted.

These correlation routines are used by the optional wall-impact coupling below.
Child and film-patch constructors are also exposed separately. Lateral wall-film
transport remains open. Inputs are checked for finite/bounded
values. Negative mass fractions produced by the reference fits return failure
and zero mass outputs; they are not clipped into a different distribution.
In particular, thick-film azimuth fits can produce negative fractions even
within the reference angle range. Velocity outputs preserve signed correlation
values; classification and child generation must interpret them together.

See [pinned source](../references/spray_wall_impact.json) and
[focused C++ comparison](validation/wall-impact-kernel-20260916.md).


### Splash children and deposited film patches

`construct_splash_children` constructs four descendants for a supplied splashed
fraction. Their droplet masses follow the four azimuthal fractions; multiplicity,
temperature and material identity are inherited. IDs start from a supplied unused
ID and `parent_id` retains the incident parcel ID. Deformation and RT timers reset.
Velocities use the reference tangent basis and splash components; children start
one child diameter outward from the impact point. A coordinate-axis fallback
handles the reference's degenerate tangent case. Negative outgoing normal speed
or invalid mass distributions reject without replacing an existing child array.

`construct_wall_film` constructs a stationary cylindrical patch from a supplied
deposited fraction, incident Weber/Reynolds numbers and contact angle in radians.
It retains liquid mass, temperature, material ID, location, normal, film diameter
and thickness. Geometry uses the pinned modified Pasandideh-Fard expression for
the coalesced parcel population. These are individual deposited patches. The
optional thermal/evaporation update below evolves each patch; merging and lateral
film motion remain open.

Both routines return a signed five-component wall budget: mass, three momentum
components and total energy, positive toward the wall. Mass remains entirely in
children or the film (wall-budget mass is zero). Momentum and energy are the
incident selected fraction's inventory minus the constructed inventory. The
caller must remove that selected fraction from the parent exactly once. The
routines leave the input parent unchanged and commit destination state only on
success. The constructors alone do not handle gas feedback, boundary crossings,
domain ownership or coupled restart; the coupled-step and application layers
below supply those connections for regular stationary walls.
See [state/conservation checks](validation/wall-impact-state-20260916.md).


### Complete standalone impact event

`resolve_wall_impact` accepts an incident parent, impact point, outward unit
normal, dimensionless impact parameters, contact angle and child ID range.
It classifies the event, generates children/deposition, updates the parent and
returns the signed wall budget in a single transaction. Rebound retains the
parent with the pinned velocity restitution coefficient; thermal breakup retains
the unsplashed parent with the reference mean splash velocity. Deposition and
spreading put all mass in a film patch. Dry splash divides mass between children
and deposition. Retained parents are placed one remaining-droplet diameter
outward from the supplied impact point. The caller remains responsible for
remaining flight time, domain checks and globally unused child IDs.

For wet splash, the pinned classifier leaves the deposition fraction at zero.
This event API explicitly deposits `1-ms`, where `ms` is its splashed fraction,
so incident mass closes with children and film. It does not entrain mass from
an existing film. This is an explicit closure extending the incomplete reference
branch, not a claim of exact wet-impact parity or physical validation.

Invalid distributions, inward child velocities and invalid film geometry reject
the entire event, preserving parent, child array and film output and zeroing the
wall budget. Successful calls replace the film/children output with this event's
results (empty when absent); a caller must append or accumulate them into its
persistent state. Regular gas coupling and checkpoint persistence are wired into
the application below; lateral film transport remains open. See
[event checks](validation/wall-impact-event-20260916.md).


### Coupled-step wall impact API

`advance_spray_coupling_3d` and `advance_spray_reactive_3d` accept optional
`wall_options` (`spray_wall_options`) and an allocated `wall_films` array.
Enable individual faces with `enabled(6)`, provide wall temperatures (K), a
contact angle (radians), and a maximum film-patch count. Enabled faces must be
slip, no-slip or isothermal gas walls. Coupled impact currently requires stationary
walls. Surface tension and viscosity come from the selected per-material breakup
properties, or global spray options even when aerodynamic breakup is disabled.

The first crossed plane is identified from the parcel substep's position segment.
Impact is applied at that substep's end; retained parents and children are placed
inside the domain and advanced for the remaining coupling interval. This shares
the current first-order particle substep approximation; exact hit-time integration
and simultaneous corner-impact resolution are not claimed. New descendants count
against existing parcel/substep limits. Invalid impact, storage overflow or later
chemistry/flow failure rolls back gas, all parcels, film patches and budgets.

Boiling temperature is inverted from the selected saturation model at sampled
gas pressure: Clausius with reference latent heat, or Antoine strictly within its
fit range. Kv uses the incident normal speed, current liquid density, viscosity
and diameter; Kcrit uses gas pressure converted from Pa to bar. Existing stationary
patches covering the impact point contribute their summed heights to the wet-film
ratio. Patches are appended individually. Their optional thermal/evaporation
update is described below; they are not merged.

`spray_system_totals(...,wall_films=...)` includes film mass, material-resolved
species inventory and integrated liquid energy. The boundary ledger includes
signed wall impulse/energy exchange; deposited liquid remains in system mass.
This capability is also exposed by the application inputs below. See
[focused coupling checks](validation/wall-impact-coupling-20260916.md).


### Wall-impact application and restart

Set `wall_impact_faces(6)` in x-low/x-high/y-low/y-high/z-low/z-high order to
enable Ahamed impacts on chosen stationary wall faces. Specify positive
`surface_tension` and `liquid_viscosity` (or per-port overrides), even with
`breakup_model='none'`. The existing `wall_temperature(6)` supplies impact wall
temperature, including on adiabatic gas walls. `wall_contact_angle` defaults to
pi/3 radians; `maximum_wall_films` defaults to one million patches. Aerodynamic
breakup can remain disabled while wall impacts generate descendants.

Impact-enabled runs write `.films.csv` with material ID, location, normal, mass,
temperature, diameter and height. History adds `film_mass`; `liquid_mass` remains
airborne parcel mass and `system_mass` includes gas, parcels and films. Parcel
CSV adds `parent_id` even without aerodynamic breakup or computational splitting.

The opt-in `_FILM_V1` checkpoint suffix stores the film count, all film records
and parcel lineage. Restart validates material identity, physical bounds,
attachment to enabled walls, mass/volume consistency and aggregate inventories.
Changes to impact settings reject restart. Film-free formats retain their prior
layout. The writer rejects nonempty supplied films unless film serialization is
explicitly enabled, preventing accidental loss of deposition state.

`cases/spray_3d/wall_splash.nml` runs a 300 m/s water parcel into the low-z wall,
producing four descendants and a static deposited patch. Liquid properties are
illustrative and gas chemistry/transport and evaporation are disabled to isolate
impact bookkeeping. See [CLI/restart evidence](validation/wall-impact-application-20260916.md).


### Stationary film heating and evaporation

Set `wall_film_evolution=.true.` and supply `wall_liquid_conductivity` in
W/(m K) for the default constant model. For temperature dependence, set
`wall_conductivity_model='polynomial'`, `wall_conductivity_coefficients(4)` and
`wall_conductivity_temperature_range(2)`. The pinned PelePhysics cubic is
`lambda(T)=a+b*T+c*T**2+d*T**3`, with T in kelvin and coefficients in SI units.
The fit must stay positive and at most 1e4 W/(m K) throughout its stated interval;
endpoints and interior extrema are checked. Film updates outside the interval
are rejected transactionally, including candidate temperatures.

For explicit injection ports, `port_wall_conductivities(16)` overrides the scalar
(default -1 inherits it). `port_wall_conductivity_models(16)` defaults to blank,
which inherits the complete global model. An explicit `constant` or `polynomial`
uses that port's `port_wall_conductivity_coefficients(4,16)` and
`port_wall_conductivity_temperature_ranges(2,16)` (both default zero).
A scalar override of an inherited polynomial requires explicit `constant`.
The global default must remain valid even when ports override it. Unused port
settings are rejected; overrides assign distinct material IDs even if vapor
species and other liquid properties match. The conductivity is evaluated at
film temperature each substep, without mixture averaging or a boiling cap.

For example, a fit can be specified with:
```fortran
 wall_conductivity_model='polynomial',
 wall_conductivity_coefficients=.1,.001,0,0,
 wall_conductivity_temperature_range=300,600,
```
These are illustrative coefficients, not a calibrated liquid dataset. See
[conductivity and restart checks](validation/wall-film-conductivity-20260916.md).

`minimum_film_height` defaults to 1e-7 m (the reference 1e-5 cm cutoff). These
inputs enter checkpoint identity. Existing `heat_enabled` controls gas/wall
heating, and `evaporation_enabled` controls film mass loss as well as droplets.
The default remains static deposition when film evolution is disabled.

For film area A and thickness h, the update uses gas conductance lambda_g/h,
wall heat `2*A*lambda_l*(Twall-Tfilm)/h`, and the pinned Chilton-Colburn mass
transfer correlation with the existing binary-carrier diffusion coefficient.
The gas skin uses the reference coefficient 1/2 for both temperature and vapor
mass-fraction interpolation, with a nonnegative hot-gas temperature difference.
The pinned source comment says one-third, but its actual `rule` constant is 1/2.
Evolving-film checkpoints from the earlier one-third implementation are rejected
by the `wall-film-evolution-v2-half-skin` context; static-film and droplet-only
contexts are unchanged. Local properties/rates are frozen over each
bounded explicit substep. Film heat capacity and liquid energy use the selected
constant or polynomial models. The footprint stays fixed; thickness is recomputed
from remaining mass and current liquid density. An evaporating film below the
cutoff transfers its complete remainder to gas and is removed from the patch list.

Wall heat is a separate signed boundary-energy contribution. Gas receives the
discrete film energy loss plus heat supplied by the wall, and vapor enters the
selected species with zero injected momentum for stationary film. This preserves
the established formation-energy convention; it is not a byte-for-byte port of
the reference energy-source bookkeeping. Gas admissibility and film property
bounds are checked before a substep commits. Newly deposited films only evolve
for their remaining coupling interval; late failure rolls back the whole step.

`cases/spray_3d/wall_film_evaporation.nml` demonstrates deposition followed by
film heating/evaporation. This remains a stationary, fixed-footprint approximation,
without lateral flow, moving-wall coupling, entrainment or calibrated film
qualification. See [focused evidence](validation/wall-film-evolution-20260916.md).


### PelePhysics drag correlation

`drag_model='pelephysics'` selects the pinned `Drag.H` correlation:
Cd is 24/Re for 0 < Re <= 1, 24/Re*(1+Re**(2/3)/6) for
1 < Re <= 1000, and 0.424 above 1000. At zero relative velocity its force
and rate are zero. With `breakup_model='etab'`, Cd is multiplied by
`1+2.6432*deformation`, using the parcel's current deformation. A negative
multiplier rejects the update transactionally. Other breakup models do not
apply this multiplier.

The default `drag_model='schiller-naumann'` retains the preceding correlation.
The selected model enters checkpoint identity, so switching drag models on
restart is rejected. The gas-skin Reynolds number, frozen exponential velocity
update and conservative mass/momentum/energy exchange are shared by both
choices. This ports the reference force law, not its complete time integration
or all spray submodels. See [focused evidence](validation/spray-reference-drag-20260916.md).


### Particle-number-weighted size distribution

Set `diameter_weighting='number'` with positive `diameter_shape` to interpret
the truncated Weibull CDF as a droplet-number distribution, matching the
sampling basis of the pinned PelePhysics Weibull model. Each parcel in a
single injection batch has the same multiplicity (represented droplet count)
and a different total mass. The common multiplicity is the requested batch
mass divided by the sum of its individual droplet masses. This retains the
prescribed total injection mass for any finite parcel count.

The default `diameter_weighting='mass'` retains equal parcel masses and
varying multiplicity. Both choices use deterministic midpoint quantiles and
the existing diameter bounds, rather than random reference sampling. With
shape zero all diameters are identical. `drop_diameter` is the Weibull scale;
it is not generally the arithmetic mean, despite the upstream `mean_dia` name.
The selector applies to initial, scheduled, multi-port and independently
scheduled port batches, and enters checkpoint identity. Switching weighting
on restart is rejected. This does not implement breakup-child generation.
See [focused evidence](validation/spray-number-cdf-20260916.md).


### Normal and lognormal diameters

`diameter_distribution='normal'` or `'lognormal'` enables the corresponding
pinned PelePhysics diameter law. `drop_diameter` is the untruncated arithmetic
mean and `diameter_std` is the untruncated arithmetic standard deviation, both
in meters. For lognormal samples the internal normal parameters are
`width=sqrt(log(1+(std/mean)**2))` and `location=log(mean)-width**2/2`.
Each port's diameter override supplies its own mean; the global standard
deviation and distribution selector provide defaults for the per-port settings below.

The existing `diameter_min` and `diameter_max` truncate and renormalize the
CDF; midpoint quantiles are inverted by bounded bisection. Positive-tail
normal probabilities use the complementary CDF to avoid subtraction from one.
Intervals with numerically unresolved probability are rejected. Positive
standard deviations must lie in [1e-15, 1] m; zero produces the mean diameter
and requires that it lie within the configured bounds. Truncation generally
changes the actual mean and standard deviation of the injected population.

Choose `diameter_weighting='number'` for upstream number-distribution semantics,
or `'mass'` for a mass-CDF interpretation. Both normalize to the requested
injected mass. `diameter_shape` must remain zero for these two laws.
The default `'weibull'` retains the previous shape-zero monodisperse or
positive-shape Weibull behavior and requires `diameter_std=0`.
Selector and standard deviation enter restart identity. See
[normal/lognormal checks](validation/spray-normal-sizes-20260916.md).


### ChiSquared table

`diameter_distribution='chi-squared'` selects the pinned 100-bin `ChiSquared`
distribution. `drop_diameter` (or each port's diameter override) maps to the
upstream nominal `d32` parameter. At bin j=1..100, its cumulative probability
is proportional to `1-exp(-xi)*(1+xi+xi**2/2+xi**3/6)` with `xi=.12*j`;
its diameter is the left edge `.12*(j-1)*d32/3`. This reproduces the source's
discrete table rather than substituting a continuous inverse distribution.
The nominal input does not guarantee the measured Sauter mean diameter after
binning, bounds, weighting and finite sampling.

Bins outside `diameter_min`/`diameter_max` are removed and the remaining
probabilities are normalized. The source's zero-diameter first bin is therefore
excluded by the positive minimum diameter. Empty support is rejected. Midpoint
quantiles select bins; repeated diameters are possible. Select number weighting
for the source sampling basis, or mass weighting for the alternative mass CDF.
`diameter_shape` and `diameter_std` must both be zero. Existing multi-port,
scheduling and restart paths carry the selector. See
[table/restart checks](validation/spray-chi-squared-20260916.md).


### Per-port diameter distributions

Each explicit injection port can override:

- `port_diameter_distributions(16)` and `port_diameter_weightings(16)`;
  blank strings inherit the global selectors.
- `port_diameter_stds(16)`, `port_diameter_shapes(16)`,
  `port_diameter_minima(16)` and `port_diameter_maxima(16)`;
  -1 inherits the corresponding global scalar.

Inheritance is independent for each field. When changing distribution, explicitly
clear an incompatible inherited parameter: for example a normal port under a
positive-shape global Weibull configuration needs `port_diameter_shapes(port)=0`.
Resolved options are validated by the same injector before any run outputs are
created. Nondefault unused-port settings are rejected. Existing port diameters
continue to set the Weibull scale, normal/lognormal mean or nominal ChiSquared d32.

Initial, shared-schedule and independent-schedule batches use these resolved
settings. Each batch conserves its assigned mass fraction or scheduled dose.
Only explicit overrides add the per-port distribution checkpoint context;
changed settings reject restart. The Fortran multi-port API accepts an optional
`spray_size_options` array, one complete record per port. This changes size
statistics, without requiring distinct liquid material IDs.
See [per-port checks](validation/spray-port-sizes-20260916.md).


### Finite circular nozzle outlets

`nozzle_radius` defaults to zero (point injection). A positive value in meters
places injection positions in the circular outlet plane perpendicular to the
specified injection axis. `port_nozzle_radii(16)` defaults to -1, inheriting the
global radius; zero explicitly selects a point outlet for that port.

Following pinned `SprayInjection.H`, the radius fraction is sampled linearly:
`r/R=(i-.5)/N`. This is uniform in radius, not uniform in disk area. Azimuth uses
the existing deterministic ID-based sequence and the same transverse basis as
the velocity cone. Each injected batch retains its previous mass and velocity
rules. The nozzle is a parcel source; it does not carve a gas inflow opening
or replace gas boundary conditions.

The application checks the entire disk's Cartesian extents before running:
`extent_j=R*sqrt(1-axis_unit_j**2)`. The lower edge may touch zero; the upper
edge must remain strictly below the domain extent. This applies even with
periodic gas boundaries. Negative/nonfinite radii and unused port overrides
are rejected. Initial, shared and independently scheduled injection all use
the resolved radii. Positive/overridden radii enter restart identity; changed
geometry rejects restart, while default point contexts remain unchanged.

The Fortran injectors accept optional `nozzle_radius` and, for multi-port use,
`port_nozzle_radii` (complete values, no inheritance sentinel). This ports
finite outlet geometry and radial sampling. Hollow and swirl options are below;
time-varying geometry and reference random trajectories remain open. See
[geometry/restart checks](validation/spray-nozzle-20260916.md).


### Solid/hollow nozzle profiles and swirl

`nozzle_profile='legacy'` retains the original equal-solid-angle velocity cone.
`'solid'` selects the source's linear relation between radius fraction q and
velocity inclination: r=R*q and theta=`cone_angle`*q. `'hollow'` places all
parcels at r=R and uses theta=`cone_angle`+`hollow_spread`*(q-.5). Both new
profiles require positive nozzle radius. Here q=(i-.5)/N, a deterministic
midpoint replacement for the source's random samples. Hollow inclination
spread defaults to zero and is allowed in [0,pi] radians; it can produce
negative or backward inclinations for sufficiently wide spreads, as in the
source. Non-hollow profiles require zero `hollow_spread`.

`swirl_angle` rotates velocity azimuth relative to the outlet position azimuth
without changing speed or axial velocity. It defaults to zero and must lie in
[-pi,pi] radians. All these application angles are radians; the upstream jet
input uses degrees. `cone_angle` remains a half-angle.

Per-port overrides are `port_nozzle_profiles(16)` (blank inherits),
`port_swirl_angles(16)` (-999 inherits), and `port_hollow_spreads(16)`
(-1 inherits). Inheritance is independent: a solid override under a hollow
global profile must explicitly clear any inherited spread. The selectors and
angles enter restart identity. Existing default contexts are unchanged.
The Fortran injectors accept optional `spray_nozzle_model` records, globally
or as a complete array per port. See
[profile/swirl checks](validation/spray-nozzle-profile-20260916.md).


### Time-dependent rate of injection

Set `injection_profile_points` to 2..64, with matching
`injection_profile_times(64)` in seconds and `injection_profile_rates(64)` in
kg/s. Times must strictly increase, cover the entire injection window and lie
in [0,1000] s; rates must lie in [0,1] kg/s. Entries after the declared count
must be zero. `mass_flow_rate` must be zero when the table is supplied.

The injector linearly interpolates between knots and integrates each scheduled
window exactly by trapezoids, including partial knot intervals. A batch is
still deposited at the beginning of its injection window, as in the existing
schedule; this is not continuous-in-time particle creation. Zero-dose windows
advance the schedule without generating parcels. Shared ports split each dose
by their configured fractions and retain their own liquid/size/nozzle settings.
The global table uses the shared schedule. Independent-port tables are configured
separately below; a global table cannot be combined with independent mode.

Table data enter checkpoint identity. Restart reconstructs both accumulated
mass and the number of nonzero batches, including pauses, before accepting
particle/event counts. Changing a knot or rate rejects restart. A zero point
count preserves the previous constant-rate behavior. See
[rate-profile checks](validation/spray-injection-profile-20260916.md).


### Independent port rate histories

With `independent_port_schedules=.true.`, each port may use
`port_profile_points(16)` and corresponding columns of
`port_profile_times(64,16)` / `port_profile_rates(64,16)`.
A zero point count retains that port's constant `port_mass_flow_rates` value.
A positive count (2..64) requires the corresponding constant rate to be zero.
Table and constant-rate ports may coexist.

Each table must cover its own `port_injection_starts`/`port_injection_ends`
window, with increasing times and nonnegative rates in the same ranges as the
global table. Unused points and ports must have zero table entries. Each event
uses the exact integral over its own port interval; zero-dose events advance
the schedule without creating parcels. Simultaneous events retain port-index
ordering. The global table and these independent tables are mutually exclusive.

Tables enter restart identity. Restart replays the merged port event order,
reconstructs individual window doses and checks aggregate injected mass and
nonzero-batch parcel counts. See
[independent profile checks](validation/spray-port-rate-profiles-20260916.md).


### Velocity from injection rate

`velocity_from_rate=.true.` derives scheduled injection speed from
`U=mdot/(rho*pi*nozzle_radius**2*discharge_coefficient)`, following the pinned
`SprayJet::readROI` relation. `discharge_coefficient` defaults to 1 and must be
in [1e-6,1]; `port_discharge_coefficients(16)` defaults to -1, inheriting the
global value. Each active configuration requires a nozzle radius of at least
1e-12 m. Density is evaluated from that port's liquid model at its injection
temperature. Derived speed is checked against the injector's existing limit.

Here mdot is the integrated dose divided by the batch window duration. Shared
ports use their assigned fractions; independent ports use their own constant
or tabulated rates and intervals. The flow-derived speed replaces prescribed
`injection_speed`/port speeds in this mode; cone, hollow and swirl settings still
set its direction. The application probes peak-rate speeds before running.

This mode requires `initial_liquid_mass=0` because an initial inventory has no
associated injection window. Coefficients are otherwise unused and cannot be
set away from their defaults when the mode is disabled. Mode and coefficients
enter restart identity. The default prescribed-speed behavior retains its
context. Each batch still has one speed and is emitted at the window start;
continuous momentum-flux equivalence is not claimed. See
[flow-derived velocity checks](validation/spray-rate-velocity-20260916.md).

## Cartesian particle partitioning

`spray_partition_mod` supplies deterministic ownership and send-buffer ordering for
the forthcoming particle MPI path. Each Cartesian axis assigns remainder cells to
the lowest coordinates. Rank numbering is x-fastest. Physical coordinates must be
inside the half-open box; periodic coordinates wrap for ownership without changing
the stored particle position. Partitions cannot exceed cell counts or one million
ranks in total.

`partition_spray_parcels` accepts globally increasing particle IDs and groups records
by destination rank, preserving ID order within each group. Counts and zero-based
offsets identify contiguous send segments. The complete parcel record and optional
per-drop component matrix move together. Invalid input leaves every caller array
unchanged. MPI communication and receive-side ID ordering are available below; distributed
gas feedback remains open. See [partition evidence](validation/spray-partition-20260916.md).

### MPI particle migration

With `PELEF_ENABLE_MPI=ON`, `migrate_spray_parcels` uses `mpi_f08` all-to-all
communication to move records to Cartesian owner ranks. IDs and lineage use integer
64-bit buffers; the 14 real parcel fields and optional component masses use separate
double-precision buffers. Received records are sorted by ID before committing. All
ranks must supply identical geometry, decomposition and component-column count.
The caller supplies globally unique IDs; duplicates arriving at one rank reject.

Invalid local input, mismatched configuration, receive-capacity overflow or invalid
received records reject collectively before any caller arrays are replaced. This
guarantee concerns application validation, not failed processes or MPI transport
errors. Empty sources and destinations are supported. This migration API does not
yet distribute gas coupling, film patches, injection IDs or application checkpoints.
See [MPI evidence](validation/spray-mpi-migration-20260916.md).

### Collective particle ID reservation

`reserve_spray_ids` reserves a contiguous interval for each rank from a shared
64-bit `last_id`. Rank order determines interval order. Each rank may request zero
to one million IDs; a zero request returns `first_id=0`. On success every rank
receives the same updated high-water mark. Negative requests, mismatched existing
marks and overflow reject collectively without changing either output value.

The caller must preserve the high-water mark in application restart state and
reserve before constructing injected or split particles. The API does not infer
that mark from locally owned records. Assignment is deterministic for a given rank
layout and request vector; changed-rank particle identity parity is not claimed.
Injection/breakup driver wiring remains open. See
[ID evidence](validation/spray-mpi-ids-20260916.md).

### Collective computational splitting

`split_mpi_spray_parcels` applies the existing multiplicity cap on local candidates,
reserves globally disjoint child IDs and commits particle/component state and the
shared ID mark together. Children inherit per-drop component masses and parent IDs;
KH initial count and shed volume are partitioned with multiplicity. If any rank
cannot split within capacity, no rank commits or consumes IDs. Input IDs must be
locally increasing and no larger than the shared mark. The component table is a
caller-owned common convention; column counts must agree across ranks.

This is computational splitting of represented drop counts, not aerodynamic
breakup. MPI injection and distributed gas/film coupling remain to be connected.
See [split evidence](validation/spray-mpi-split-20260916.md).

### Collective cone injection

`inject_mpi_cone_parcels` combines global ID reservation, the existing cone injector
and a collective commit. Each rank supplies its own injected mass, count, position,
velocity and size/nozzle options. Zero-count ranks supply zero mass. Optional common
component tables and rank-local fractions create mixed droplets using density and
energy evaluated at injection temperature. Component columns append with particles.
Returned `injected_totals` are local mass, momentum, energy and gas-species amounts;
the caller sums or records them as appropriate.

A failed local geometry/property check or capacity check preserves every rank's
particles, component matrix and shared ID mark, returning zero injection budgets.
Successful injection appends increasing global IDs but does not automatically move
particles to spatial owners; call migration after injection. Distributed scheduling
and the MPI CFD application remain open. See
[injection evidence](validation/spray-mpi-injection-20260916.md).

### Distributed gas-field sampling

`sample_mpi_spray_field` takes each rank's Cartesian cell block and local particle
positions. It sends the eight interpolation-cell indices to their owning ranks and
returns only the requested field values, then evaluates the ordinary particle
stencil locally. No root or replicated whole-grid field is required.

All ranks must agree on global cells, spacing, lower corner, Cartesian parts,
boundary faces and field width. `spray_cartesian_block` supplies zero-based first
cells and local extents. Empty requesting ranks are supported. This bounded API
accepts up to 125000 local positions, 256 fields and one million incoming cell
requests per rank. Invalid configuration, positions, field values or received
indices reject collectively while preserving previous sample arrays. These are
interpolated conserved fields; thermodynamic recovery remains the caller's job.
Distributed feedback deposition and time stepping remain open. See
[sampling evidence](validation/spray-mpi-sampling-20260916.md).

### Distributed feedback deposition

`deposit_mpi_spray_field` distributes local particle inventory increments through
the same eight-cell stencil. It routes indices and weighted contributions to cell
owners, divides by cell volume and adds to local field candidates. It needs no root
or whole-grid copy. Sources are accumulated in sending-rank order, then local
particle/stencil order; changed-rank floating-point parity is not implied.

The geometry, field-width and request limits match sampling. Empty source ranks
are supported. Invalid input, request routing or floating-point overflow rejects
collectively and preserves every local field. This is a generic conservative
field operation: gas positivity, species closure and temperature recovery must be
checked by the future coupled driver before committing its particle candidates.
See [deposition evidence](validation/spray-mpi-deposition-20260916.md).

### Thermodynamically checked distributed feedback

`deposit_mpi_spray_gas` applies inventory increments to distributed conserved gas
blocks and recovers temperatures with the existing reactive thermodynamics. It
validates the initial state, deposits into a private candidate, and accepts both
state and temperature only when every owner recovers successfully. Invalid density,
species or thermal energy on one owner preserves both arrays on all ranks. The
ordered species table must be identical on every rank. Geometry and request limits
are inherited from field deposition.

The caller still owns particle candidates, time-step retries and the combined
particle/gas commit. This API does not itself advance particles or chemistry.
See [gas feedback evidence](validation/spray-mpi-gas-20260916.md).
