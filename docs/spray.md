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
spatially isothermal, with constant or temperature-dependent polynomial density. `vapor_name` must identify a gas
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
multicomponent liquid properties or temperature-dependent heat capacity.

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
creating liquid mass. Liquid heat capacity and the existing formation-energy
reference remain unchanged; liquid p/rho and surface-energy corrections are
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
those fractions do not multiply the independent scheduled rates. Material,
velocity magnitude, temperature and size settings remain shared.

Restart preserves the aggregate event count and reconstructs the per-port
counts by replaying the deterministic schedule. It checks injected mass,
parent counts, future events and missed events before continuing. Schedule
arrays are included in checkpoint identity. Replay is bounded by the input's
one-million-event limit; large restart performance is not qualified.
See `cases/spray_3d/water_independent_ports.nml` and
[focused evidence](validation/independent-port-schedules-20260915.md).
