# Completion roadmap

## Definition and frozen reference

A selected PeleC responsibility is complete only when it has an independent
Fortran implementation and automated numerical evidence, or is explicitly
excluded with a reason and a usable alternative. A scaffold, test count or
clean compilation is not sufficient.

The reference is PeleC `development`
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`; recursive dependencies are pinned in
[the manifest](../references/pelec_baseline.json). Changing that baseline must
identify affected numerical signatures and acceptance criteria.

## Consolidation versus solver completion

The source-provenance repair, documentation ownership, CMake module split,
archived PR11 supersession review, and hosted serial/MPI/CVODE/clean-install
workflows are implemented. [PR #102](https://github.com/GTT10/Fortranslate/pull/102)
records the actual integration status and completed run identities; see
[current status](current_status.md) and [the layout evidence](validation/build-layout-20260909.md).
A targeted NASA7 failure-status fix and its red/green tests are recorded
[separately](validation/nasa7-conversion-status-20260909.md).

Do not redo this housekeeping as a substitute for the remaining physical
capabilities. Conversely, do not label the following capabilities implemented
merely because the repository has been consolidated.

## 0.246.0 functional increment

The opt-in regular 3D path now connects coexisting single-component liquid populations selected per port, conservative
mass/momentum/total-energy exchange, drag/heating/evaporation, timed cone
injection from point or finite circular solid/hollow outlets with swirl with up to 16 independently configured ports, shared or independent constant/tabulated-rate schedules with prescribed or discharge-based speeds and optional per-port mass- or number-CDF Weibull/normal/lognormal/ChiSquared diameters, deviatoric Smagorinsky and coupled checkpoint/restart to the existing
gas solver. Cell-local CVODE is connected to this CFD path, with a pinned
FFCM1_Red liquid-methanol burn fixture and independent Cantera gas-only oracles.
Chemical labels, single-interval NASA7 and classified Lindemann/SRI/Tsang, chemically activated, PLOG and Chebyshev ingestion are
supported. Selected mechanism capacity is now 128 species; a 53-species GRI-Mech
0D and coupled-spray case has focused evidence. Large-mechanism MPI/scaling and
production diesel chemistry qualification remain open. Lu n-dodecane (53 species,
268 reactions) now imports with extended Troe coefficients and completes a
lean constant-volume ignition comparison; Newton trial-state handling preserves
accepted-state checks. See [focused evidence](validation/dodecane-20260915.md). See [scope and commands](spray.md) and [acceptance gates](validation/0.246.0.md).
These are bounded implemented capabilities, not completion of diesel chemistry,
multicomponent liquids, general dispersion/collision, characteristic inlet models and near-wall turbulence qualification or MPI/AMR particles.
Regular wall reflection and outflow now include gas/particle surface budgets
and restart persistence; isothermal/tangentially moving gas walls include molecular
and SGS heat/work budgets. The Ahamed impact-regime, splash-velocity, mass-distribution, four-child construction and deposited-film geometry kernels are connected by a transactional impact event with explicit wall impulse/energy budgets; regular coupled-step boundary handling now retains deposition patches, advances splash children and evolves stationary film temperature/mass with per-material constant or cubic conductivity; the CLI and coupled checkpoint retain these states, while lateral film transport and physical qualification remain open. ETAB and KH-RT representative-drop
size/deformation, bounded computational parcel splitting and restart are implemented;
KH shed-mass accumulation and conservative four-child transverse release now
connect to CFD and restart, including computational splitting with partitioned
KH inventories and stationary wall impact processed before KH release on a
crossing substep. General turbulent dispersion, wet-impact correlation limits
and physical qualification remain open. Uniform gas/particle acceleration
now retains external impulse/work separately through coupled restart. Single-component
Antoine saturation-pressure fits are connected with explicit temperature bounds
and restart identity. Polynomial liquid density now drives injection, drag,
evaporation and breakup consistently; inverse-temperature liquid viscosity now
drives ETAB/KH-RT. The pinned piecewise drag law and ETAB deformation multiplier are selectable alongside Schiller-Naumann drag, with restart identity. An optional near-boiling vapor-state temperature cap follows
the pinned reference convention. Port material identity now routes evaporation to each liquid's gas species and persists through restart; density/pressure fits, vapor-state caps and breakup viscosity/surface tension can be set independently per port. Polynomial liquid heat capacity now supplies the thermal update and integrated liquid energy, including per-port fits. Mixed-liquid density, heat capacity, stored energy and fixed-composition heating
now have a shared thermodynamic API. Raoult surface partitioning and supplied
component-mass evaporation exchange, diffusion-weighted rates and a conservative
frozen-gas transport update are implemented. Regular-grid feedback, persistent
component masses, boundary budgets, representative breakup/splitting and composition
checkpoint APIs are connected. Combined hydro/transport/LES/chemistry stepping
now retains mixture composition through accepted updates and rollback. Physical mixture boiling rates remain open.
The mixture impact event kernel now partitions component masses into retained
parent, splash and film with wall budgets. Supplied film evaporation/heat exchange
now updates composition, temperature and height. Mixed film heat/evaporation rates
and a frozen-rate update are connected to grid feedback. Planar
crossing dispatch now computes mixed impact parameters and retains component budgets.
The checkpoint API retains mixed-film component masses. Stationary mixed-film grid
feedback now updates gas and film with wall budgets. Particle crossing, splash and
remaining-interval film evolution now share one transaction. Reactive stepping
retains both component matrices with full rollback. CLI wall inputs, component-film
output and restart are now connected. Shared-composition CLI inputs, initial/scheduled injection,
component output and restart are now connected. Mixed KH shed children preserve
composition and conservative transverse energy through coupling, splitting and restart.
Component density/heat-capacity fits, Antoine pressure and Clausius surface caps
are connected to CLI inputs and restart identity. Per-port compositions now feed
initial/scheduled injection, component budgets and density-based injection speeds.

## Ordered implementation work

| Order | Work unit | Measurable exit gate |
| --- | --- | --- |
| 1 | Validate existing physics before expanding it | Pin external cases, states, observables and tolerances; compare H2/O2 trajectories/rates with Cantera and flow fields with frozen PeleC; characterize convergence and conservation |
| 2 | Production chemistry path | Specify the target fuel mechanism and unsupported forms; validate ingestion/0D ignition; qualify stiff CFD integration with rollback and conservation before detailed-fuel claims |
| 3 | General 3D flow and mesh | Physical boundaries, then dynamic/multipatch/multilevel AMR and general EB in separate increments; require conservative regrid/reflux, changed-rank restart and convergence |
| 4 | LES and Lagrangian spray | Validate closures and particle submodels, then two-way mass/momentum/energy coupling; qualify evaporation and breakup before penetration/ignition comparisons |
| 5 | Practical release qualification | Reproducible clean install, documented supported inputs, distributed I/O/scaling measurements, physical benchmarks and owner-resolved licensing |

These are dependencies, not permission to bundle several physical models into
one PR. Existing regular 3D and static two-level AMR foundations must not be
reimplemented. Their restrictions remain explicit in current status.

Regular 3D hydro/chemistry/molecular transport now accepts per-face prescribed reservoir inflow, outflow, slip-wall and adiabatic/isothermal no-slip conditions with prescribed tangential wall velocity
with integrated boundary-flux accounting; see the [case guide](../cases/reactive_boundaries_3d/README.md).
Characteristic inflow and moving geometry, general EB and AMR-particle boundary coupling
remain part of the completion work.

Serial/MPI 3D temperature-tagged patch movement/resizing is implemented with conservative
remapping, overlap reuse and geometry-aware changed-rank restart, including
physical-boundary tagging and boundary-touching patches. Remapping now runs
within fixed coarse slabs without full-level root aggregation. General
patches/levels, measured workload/scaling qualification, distributed checkpoint I/O
and EB coupling remain open. Two-level MPI AMR CSV output and final diagnostics now
operate on local slabs, as does analytic coarse/fine initialization; see the [AMR case guide](../cases/amr_reactive_3d/README.md).

Particle Cartesian MPI migration now preserves parcel/component payloads and
restores local ID ordering, with collective input/capacity rejection. Distributed
gas feedback, film ownership and application restart remain. A collective ID
reservation, collective computational splitting and cone injection are implemented;
distributed scheduling and aerodynamic-breakup driver wiring remain. Owner-routed
particle field sampling and feedback deposition are implemented; coupled
thermodynamic stepping remains.

## Contract for every increment

Before modifying code, specify the failure/capability, supported dimensions,
physics/parallel mode and a measurable exit test. Include a focused unit test,
an application/conservation or external-parity test, and a precise scope
statement. Separate numerical changes from structural refactoring.

For comparisons, record reference revision, mechanism/source hash, initial
state, observable, tolerance, comparison time and mesh/time-step convergence.
Do not loosen tolerances or regenerate golden outputs merely to obtain green
CI. Failed, interrupted, missing and completed checks are distinct states.
New work starts from the integrated main; a deliberate PR dependency must
state its base and merge order rather than recreating an accidental stack.

## Historical detail

The old milestone-by-milestone roadmap remains at the
[immutable checkpoint](https://github.com/GTT10/Fortranslate/blob/509d07929d52d9e570a6af1f6ec7e5836bd4fe5d/docs/completion_roadmap.md).
Detailed rationale remains in [design decisions](design_decisions/), the
[historical implementation inventory](implementation_status.md) and dated
[validation records](validation/). [Issue #101](https://github.com/GTT10/Fortranslate/issues/101)
tracks completion and actual tooling limitations. The unexecuted ancestor-only
branch-deletion audit is an administrative remainder, not a numerical feature.
