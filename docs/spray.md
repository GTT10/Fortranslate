# Experimental periodic 3D spray combustion (0.246.0)

This is a working opt-in Eulerian--Lagrangian path, not production diesel-spray
qualification. It couples the existing regular periodic 3D gas solver to
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
Particle locations wrap periodically. Gas interpolation and source deposition
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
boiling/real-fluid handling, multicomponent liquids, breakup/collision, calibrated size distributions
and turbulent dispersion, physical injector/wall boundaries, particle
migration with MPI/AMR/EB, scalable output/checkpointing, detailed diesel-fuel
mechanisms and performance tests. Selected chemistry retains its 2--32-species
limit and currently supported reaction forms; FFCM1_Red does not represent
HXN/HMN/AMN. No claim of measured ignition, penetration or liquid length is made.
