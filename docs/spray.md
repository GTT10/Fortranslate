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

SI units are used throughout. The liquid is dilute, spherical, single-component,
constant-density and spatially isothermal. `vapor_name` must identify a gas
species in the selected mechanism. Liquid property numbers in the examples
are illustrative constants, **not** a fitted methanol/water property database.

The reference energy is `e_l(Tref) = h_v(Tref) - Lref`, using the selected
mechanism's NASA7 gas enthalpy. Thereafter `e_l(T) = e_l(Tref)+cp_l*(T-Tref)` and
the latent heat in the droplet energy equation is `h_v(T)-e_l(T)`. The
incompressible liquid pressure-work correction and surface energy are neglected.
This choice aligns chemical formation-energy references across phases. The
gas receives minus the parcel's discrete change of mass, momentum and **total**
energy; latent heat, evaporated kinetic energy and drag work are not added a
second time.

The surface saturation pressure uses a **constant-Lref** Clausius--Clapeyron
law anchored at `(Tref, pref)`. States with `psat >= pgas` are rejected rather
than extrapolated as boiling/supercritical flow. Condensation is not modeled;
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
mechanisms and performance tests. Selected chemistry retains its 2--32-species
limit and currently supported reaction forms; FFCM1_Red does not represent
HXN/HMN/AMN. No claim of measured ignition, penetration or liquid length is made.


## Regular-grid wall and outflow boundaries

`boundary_faces(6)` in `&spray` uses x-low, x-high, y-low, y-high, z-low,
z-high order. The default remains periodic. Supported alternatives are
`inflow`, `outflow`, `slip_wall`, and stationary adiabatic `no_slip_wall` for gas; parcels
reflect specularly at either wall type. Opposing periodic faces must be paired.
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
moving/isothermal spray walls, AMR/EB/MPI parcels and physical wall-impact models
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


## KH-RT representative-drop breakup

Set `breakup_model='khrt'` with positive `surface_tension` and
`liquid_viscosity`. `khrt.nml` uses the pinned PelePhysics defaults
`khrt_b0=.61`, `khrt_b1=7`, `khrt_c3=1`. KH waves relax the representative
radius; RT instability uses a persistent waiting time, initialized to -1
and enabled by the first KH event. The Reynolds number uses film density and
viscosity; the wave equations use bulk gas density and relative velocity.
Mass reduction increases multiplicity at fixed represented mass, velocity and
temperature. The minimum diameter and transactional failure rules still apply.

This implements the upstream `make_new_drops=false` branch. Shed-mass children
of a distinct radius, transverse velocity dispersion and primary atomization
remain open. Computational splitting below partitions identical representative
droplets. Liquid properties are prescribed constants; surface energy is omitted.
Zero RT acceleration resets an active timer to avoid singular expressions.
Exact upstream timestep parity and physical validation have not been established.

Parcel CSVs include `rt_time`. KH-RT checkpoints use V7/V8 (periodic/physical)
or V9/V10 with computational splitting, retaining the RT timer and binding all
three coefficients into the restart context. Existing V1--V6 formats are unchanged.
See [source provenance](../references/spray_khrt.json) and
[focused evidence](validation/khrt-20260915.md).


### Computational parcel splitting

For ETAB or KH-RT runs, set `maximum_multiplicity` to a positive upper bound on represented
droplets per computational parcel. Zero (the default) disables this partition.
`maximum_parcels` bounds allocated particle storage, up to 1,000,000.
After a coupling interval, overweight live parcels are partitioned into equal
weights, with the final child carrying the rounding remainder. Droplet mass,
position, velocity, temperature, deformation and RT timer are copied; aggregate mass,
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
See [splitting/restart evidence](validation/spray-splitting-20260915.md).


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
