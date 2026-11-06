# Reacting Taylor–Green AMR reference control

This case uses the unchanged pinned PeleC `TGReact` problem with LiDryer,
MOL with two iterations and CVODE chemistry. The periodic 1 cm cube has a
4³ root and a ratio-two fine patch covering its lower x half. It advances
three 100 ns coarse steps with two fine steps per coarse step. Molecular
transport is disabled.

`profile.csv` supplies 33 points including both periodic endpoints:
`T = 2000 + 50 cos(2πx/L)` K, `Y_H2 = 0.06 + 0.002 cos(2πx/L)`,
`Y_O2 = 0.5`, the remainder N2, and pressure 101325 Pa. TG velocity amplitude
is 10 m/s. The original problem interpolates this input; the Fortran probe
imports the resulting initial cell values rather than approximating that
initializer. Input coordinates and pressure are SI; the native problem
converts them to CGS. Solver input geometry and velocity are CGS.

The imported state includes conserved gas, temperature, actual UEINT and
the initial reaction history. Native `PeleC::post_init` initializes the last
of these over the level timestep before the first MOL predictor; replacing
it with zero does not reproduce the same initial-value problem. Only time
zero is imported. Later states and reaction histories are computed by the
production Fortran AMR integrator and its CVODE backend.

Add `--initialize-history` to discard the imported reaction estimate and
compute it independently with the production Fortran initializer. This mode
still requires exact initial gas/temperature/UEINT identity and applies the
same trajectory, conservation and chemical-change acceptance limits. Initial
reaction-rate differences are recorded separately.

Build `PeleC-TGReact`, the plotfile exporter, and the Fortran
`test_amr_tree_reacting_reference` target in a tests/SUNDIALS-enabled build.
The latter uses its own SHA-checked LiDryer fixture, independent of the
user-selected application mechanism. Run:

```sh
python3 tools/check_tgreact_amr_reference.py \
  --pelec /absolute/path/to/PeleC-TGReact \
  --reference-source /absolute/path/to/pinned/PeleC \
  --probe /absolute/path/to/test_amr_tree_reacting_reference \
  --exporter /absolute/path/to/pelef_plotfile_csv \
  --work build/tgreact-amr-reference
```

Replace `--pelec` and `--reference-source` with
`--reference /absolute/path/to/captured/reference` to reuse SHA-verified
native plotfiles. `checks.json` records commands, binary and input hashes,
all four physical times, field errors, inventories and failed criteria.
The checker returns nonzero if any acceptance condition fails.

Use `--subcycling 0` for a common timestep on both levels. The default
`--subcycling 1` preserves the recursive case above. The common-clock variant
changes only the copied native input's `amr.subcycling_mode` to `None` and
the probe's subcycling flag; all field, inventory and reaction-activity limits
remain unchanged. Generate a matching native reference first: the checker
rejects a captured reference with a different timestep mode. Older manifests
without the mode field belong to the original recursive case. The selected
mode is recorded in both `checks.json` and new reference manifests.

Add `--transport --reference-source /absolute/path/to/pinned/PeleC` to
enable viscosity, thermal conduction, species diffusion and species enthalpy
transport on both sides. This imports the exact SHA- and Git-blob-checked
LiDryer Simple fits already used by the EB comparison. The native copied
input changes the four `diffuse_*` flags to one; the probe's V3 input carries
the same fits. The original V2 no-transport input remains supported. New
manifests record transport enablement and coefficient provenance, and reject
references from the other transport mode. The geometry, duration, chemistry
settings and acceptance limits remain those of the original case.

`--cells 8` or `--cells 16` changes the root mesh from the default 4 cells
per direction while preserving the physical domain and lower-x-half refined
region. Native grid inputs and the probe's V4 hierarchy use the same size;
reference manifests reject a different root mesh. The original 33-point
profile, timestep, duration and acceptance thresholds remain fixed. This
provides matched grid comparisons; a pass on a finer grid alone is not a
mesh-convergence result, and profile-sampling effects must be distinguished
from evolution errors in a convergence study.

For operator diagnosis, `--reaction-off` disables chemistry in both solvers
and requires exactly zero imported reaction history. It checks hydrogen
conservation instead of a minimum hydrogen-consumption fraction; other limits
are unchanged. It cannot be combined with `--initialize-history`.
`--mol-iterations 1` selects one MOL iteration on both sides; the default is
two. Both controls are recorded and checked when reusing references. Prefer
MOL1 for reaction-on/off isolation because the pinned nonreacting MOL2
register/conservation conflict is a separate unresolved comparison.

Acceptance is declared in `acceptance.json`: initial identity, field
agreement, nonuniform initial conditions, species closure/positivity,
mass and H/O/N conservation, numerical-energy accounting, and resolved
hydrogen consumption. The case compares the serial reacting AMR operator;
it does not qualify application initialization, regridding, MPI, EB,
transport, convergence or physical ignition.

For larger meshes, `--command-timeout 900` allows up to 900 seconds for each
native/probe/export subprocess. The default remains 120 seconds. This changes
only the wall-time limit, not physical time, integration steps or acceptance
thresholds; command records include the limit and explicitly mark timeouts.

`--time-refinement 2` or `4` divides the 100-ns coarse timestep by that factor
and multiplies the three coarse steps, preserving the 300-ns final time.
Both native and candidate runs use the selected clock. Captured references
record the factor and reject mismatched reuse; older captures default to 1.
Original field and inventory thresholds remain unchanged.

With MPI enabled, `test_mpi_amr_tree_reacting_reference` builds from the same
probe source. Pass it as `--probe` with `--mpi-ranks 2` or `3` (and optionally
`--mpiexec PATH`). The two patches are owned by ranks 0 and 1; rank 2 is empty.
Stepping and optional history initialization use the production MPI routines.
Only output collection replicates fields for writing the serial CSV schema.
This direct reference comparison is distinct from application/restart testing.

`--x-boundary slip_wall` sets both x faces to PeleF `slip_wall` and native
`SlipWall`, with x periodicity disabled and y/z periodicity preserved. The
boundary mode is recorded in reports and captured-reference manifests;
mismatched reuse is refused. Periodic remains the default. Wall comparisons
must pass their own field/inventory checks and do not inherit periodic results.

Probe budget CSVs also record cumulative boundary and external-source integrals
for density, three momenta and total energy. The checker reconciles each with
composite inventory change and the energy correction. Signed-momentum residuals
are normalized by the largest of absolute inventory, initial/final signed
inventory and net supplied change; this avoids a zero-denominator cancellation
artifact. The existing inventory tolerance applies.

`--x-boundary no_slip_wall` instead selects stationary adiabatic `NoSlipWall`
in the native solver and `no_slip_wall` in the candidate. This reflects both
normal and tangential velocity; no moving-wall or imposed-temperature condition
is included. It uses the same field and cumulative boundary-ledger checks.

`--x-boundary outflow` selects first-order extrapolation on both x faces:
PeleF `outflow` and the pinned native input keyword `FOExtrap` (internal
boundary index 2). It imposes neither pressure nor a nonreflecting condition.
These runs additionally export species boundary/source integrals and reconcile
candidate element inventories with the supplied amounts. Native element/mass
inventories are compared with the candidate; native boundary flux is not
independently exported, so this is not a native conservation audit. H2 consumption
is adjusted using candidate boundary H2 supply for both solvers, with that basis
recorded explicitly. Closed-domain checks and their limits remain unchanged.

`--initial-velocity zero` is a separate initial-flow diagnostic. It changes
only native `prob.u0` to zero; the candidate imports that native initial state,
with exactly zero momentum required. The frozen temperature/composition profile
and native field/budget tolerances remain unchanged. The default `tg` retains
the original minimum-speed check. Captured references record the mode and
reject mismatched reuse. This removes all initial TG flow, not only the normal
boundary derivative, and does not make the piecewise-linear thermal profile
zero-gradient compatible.


The optional `--cells` values are 4, 8, 16 and 32; default4 is unchanged.
These preserve the physical domain, half-domain refinement, original33-point
profile and acceptance limits. The n32 extension permits a finer interface
resolution diagnostic; its availability does not establish reference or
spatial-convergence qualification. The fixed profile remains piecewise
interpolated, so this is not an analytic smooth-profile convergence test.
