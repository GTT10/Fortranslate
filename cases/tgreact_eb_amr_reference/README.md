# Nonuniform reacting EB-AMR reference case

This case compares the production Fortran EB tree integrator with the
unchanged pinned PeleC TGReact problem. It uses a 1 cm cube, an 8³ root,
one 8³ fine patch covering coarse indices 1:4 in all directions, a planar
wall at z = 3.3 mm, periodic x/y and adiabatic slip z boundaries. The gas has
spatially varying temperature, H2 composition and Taylor–Green velocity.
Three 100 ns root steps use recursive MOL2 and LiDryer CVODE chemistry.

Molecular viscosity, conduction and species/enthalpy diffusion are enabled.
The Fortran probe reads Simple transport fits from the exact pinned
PelePhysics header through the SHA/Git-blob/species-order checked importer;
it does not substitute the Lennard-Jones transport approximation. The input
uses the native default computable solid-body state (`eb_zero_body_state=0`).
PeleC requires EB blocking factor >= 8. NoRedist is used to isolate the
unredistributed spatial operators; this is not a small-cell stability claim.

Build `test_eb_tree_reacting_reference` in a selected-runtime build with
SUNDIALS enabled, then run:

```sh
python3 tools/check_tgreact_eb_amr_reference.py \
  --pelec <pinned-PeleC-TGReact> --reference-source <pinned-source> \
  --probe <test_eb_tree_reacting_reference> \
  --exporter <pelef_plotfile_csv> --work <artifact-directory>
```

The checker imports only the native time-zero gas/T/UEINT values. It computes
the initial reaction history independently and evolves all later fields with
the production EB reaction/hydro/transport driver. The imported history is
retained only for reporting initial-estimate differences. Covered coarse
cells and solid cells do not contribute to composite inventories; analytic
and exported volume fractions must match at every level and time.

`acceptance.json` freezes input hashes and the existing regular TGReact
field/conservation/consumption gates. A completed captured reference can be
reused with `--reference <reference-directory>`; all captured files, source
contract and input contract are checked. `--transport-off` creates a matched
intervention and is separately recorded in its capture manifest. Every
attempt retains logs, executable/input identities, physical times and
`checks.json`, including failures.

A failed comparison is an unresolved discrepancy, not authorization to
change reference algorithms or widen tolerances. This probe does not qualify
ordinary application initialization, MPI, regridding, mesh/time convergence,
ignition or experiment.

The first completed comparison fails the unchanged spatial field gates;
see the [dated record](../../docs/validation/tgreact-eb-amr-reference-20261019.md).
Conservation and close chemical consumption are not sufficient for parity.

## Paired-pressure comparison

Add `--paired-pressure` to the checker to evaluate the opt-in pressure policy
against the same unchanged captured reference and field/conservation gates.
This forwards `--paired-reacting` to the probe, enabling separate momentum,
consistent velocity and paired pressure in the production reacting and
molecular-transport MOL driver. Initial reaction history is still computed
independently. The ordinary two-argument probe keeps all three options off.
Use `--momentum-only` or `--velocity-only` as mutually exclusive checker
options to isolate the prerequisites on the same case. These map to
`--momentum-reacting` and `--velocity-reacting` on the probe.

The checker records this as an intervention with `qualification=false`, even
if every comparison passes: the policy changes the pinned spatial algorithm.
The existing `--audit-native` mode is a separate chemistry/transport-off
pressure diagnostic and cannot substitute for this reacting comparison.

## Analytic shear diagnostic

The probe's `--audit-shear` mode replaces imported fluid states with exact
cell averages of steady inviscid sinusoidal shear. The V4 input still supplies
geometry and numerical controls. `tools/check_eb_shear_reference.py` compares
all four halo policies and records conserved-field, held-energy, temperature
and pressure norms against the analytic solution. It does not use the
reacting reference as truth. See the
[case definition and artifacts](../../docs/validation/eb-shear-analytic-20261022.md).

## Bounded temporal diagnostic

`tools/check_tgreact_eb_temporal_study.py` reuses an identified completed
production capture and the separately identified native held-EOS diagnostic.
It runs 50 ns and 25 ns root steps to 300 ns, plus 100 ns root steps to
1.2 microseconds, on the same mesh. Each run independently initializes its
reaction history for its own time step. The 100 ns production capture supplies
the third resolution for volume-weighted temporal differences at 300 ns.

```sh
python3 tools/check_tgreact_eb_temporal_study.py \
  --baseline <completed-production-capture> \
  --policy-record <held-policy-artifact-directory> \
  --probe <test_eb_tree_reacting_reference> \
  --native-held <identified-held-EOS-native-binary> \
  --exporter <pelef_plotfile_csv> --work <artifact-directory>
```

The three schedules run concurrently with one thread per executable. Logs,
the checker copy, command times and per-run results remain in a unique output
directory. `status: completed` means the calculations and conservation checks
finished; inspect `diagnostic_limit_exceedances` for field and consumption
differences. Every report has `qualification: false`: using modified native
spatial EOS cannot close the unchanged pinned-reference gate. Temporal
differences on this fixed mesh also do not establish spatial convergence,
ignition or experimental validation.

## Explicit native activation basis

Pinned CEPTR and the canonical Cantera importer use slightly different calorie
gas constants when forming Arrhenius activation temperatures. The canonical
`mechanisms/lidryer.yaml` and its bundle retain the Cantera/SI interpretation.
For an isolated native-compatibility investigation, generate a separate bundle:

```sh
build/cantera-venv/bin/python tools/prepare_lidryer_native_activation.py \
  --reference-source <pinned-source> --output <new-output-directory>
```

The tool verifies the pinned YAML/header hashes, checks all 17 nonzero native
activation coefficients against the header's calorie constant, then creates
an explicit Kelvin-activation YAML source. The pinned importer regenerates its
bundle and Fortran module, including source-provenance verification. Thermo,
species, transport and every non-activation reaction parameter must remain
identical. Every activation entry must match the source-derived conversion.
The original mechanism and
reference header are untouched; an existing output directory is refused.

`activation-basis.json` records the source identity and commands. Keep this
variant's comparisons separate from canonical Cantera chemistry and from the
original pinned-PeleC gate. A held-EOS native comparison still changes the
native spatial algorithm, even when the activation basis is matched.

## Mesh diagnostic inputs

The probe accepts `PELEF_REACTING_EB_REFERENCE_V2`, followed by one integer
root-size row (`8`, `16` or `32`), then the existing V1 stepping, spacing,
plane/transport, 63 fit rows, record count and cell records. Both levels have
N-cubed cells; the fine patch covers coarse indices 1:N/2 in each direction
with ratio two. Exactly `2*N**3` unique cell records are required. The physical
domain remains N times the supplied root spacing. V1 retains its 8-cubed
layout. This is a bounded refinement family, not arbitrary patch topology.

V3 adds a `relative absolute_density_kg_m3` row immediately after the root
size. It accepts finite positive chemical integration tolerances with
relative tolerance at most `1e-10` and density absolute tolerance at most
`1e-15 kg/m3`. V1/V2 retain their historical relative tolerance `1e-10`
and mass-fraction absolute tolerance `1e-14`, which the production driver
multiplies by local density. V3 bypasses that multiplication at the reactor
callback so its explicit density tolerance reaches CVODE unchanged.

V4 adds a `0` or `1` row after the V3 reaction-tolerance row. `1` selects
PeleF's explicit pinned EOS convergence mode: a `1e-12 K` caloric Newton
stopping criterion and 5000 iterations, corresponding to the native
`CONVERGENCE` option. Physical and bounded signed-trial inverses use the same
setting. `0`, V1, V2 and V3 retain the pinned default (`1e-6 K`, 200 iterations).
The EOS policy identity includes the convergence mode, so checkpoint and MPI
policy comparisons distinguish it. Species coefficients and extrapolation
bounds do not change.

The probe reuses successful initialization reactions only when all input
bits and the step-limit option match. The cache is bounded to 128 entries,
belongs to this one immutable loaded mechanism and is not used during
evolution. A full cache falls back to a new solve; failed reactions are never
cached. It reports initial cache counts and the effective CVODE density
tolerance range. This avoids repeated identical initialization solves without
approximating states or changing production integrator behavior.

The existing checker accepts `--mesh-study <configuration.json>` with a live
`--pelec` diagnostic. The configuration declares `schema_version: 1`, `cells`,
`steps`, `dt_s`, `timeout_s`, and absolute paths for `mechanism_bundle`,
`mechanism_record`, `native_policy_identity`, and `probe_build_record`.
The mechanism record is the activation preparation tool's `activation-basis.json`;
the native identity is the held-EOS diagnostic's `identity.json`. The probe
build record binds `native_activation_binary_sha256` and
`mechanism_bundle_sha256` to the supplied executable and generated bundle.
Source/case pins and the activation-only conversion are checked before running.

An optional `reaction_tolerances` object contains `relative` and
`absolute_density_kg_m3`. Study mode then writes V3 and sets native
`ode.rtol` and `ode.atol`, converting the density absolute tolerance to
g/cm3 by multiplying by `0.001`. Typical-value scaling is disabled explicitly.
This option permits only tighter controls than the original native density
tolerance (`1e-15 kg/m3`). Native applies the same scalar absolute tolerance
to its temperature coordinate; PeleF instead recovers temperature from its
reduced-density state and energy. These formulations remain distinct.
Reports record both effective controls, and the mesh analyzer rejects mixed
chemical tolerances. The frozen field acceptance limits are unchanged.

With explicit `reaction_tolerances`, optional `fortran_eos_convergence: true`
writes V4 and selects the tighter inverse in PeleF only. The native diagnostic
binary remains at its recorded default inverse precision. Reports mark both
EOS modes, and mesh comparison rejects mixed precision controls. This is a
declared one-sided diagnostic intervention, not unchanged-PeleC qualification.

For a separate spatial-convergence diagnostic, optional `profile_points: 129`
replaces only the working copy of the 33-point table. The generator
`tools/prepare_tgreact_profile.py` verifies the frozen table against its declared
cosine temperature/H2 formulas, writes 129 points and records both hashes.
All cell centers in the 8/16/32 root/fine family then sample those analytic
formulas directly. The fixed 33-point table instead introduces interpolation
kinks whose alignment changes with mesh size. This distinction matters when
interpreting observed order. The original profile, case pins and tolerances
are preserved. The checker records the separate sampling convention and actual
profile hash; mesh comparison rejects mixtures of sampling conventions or data.
257 and 513 points are also accepted, with exact fine-center alignment required.

Study mode resizes the native root and fine box, imports the independently
sampled native initial state, checks analytic EB fractions at both levels,
and applies the unchanged field, species, element and energy-budget limits.
Success is `status: diagnostic_passed` with `qualification: false`. Study
captures are explicitly refused by ordinary pinned-reference mode. A successful
comparison on one finer mesh is not a convergence-order result; multiple
resolutions at controlled time steps and consistent volume restriction are
still required.

Compare two or three completed study directories with:

```sh
python3 tools/compare_tgreact_eb_meshes.py \
  --runs <grid8-run> <grid16-run> <grid32-run> --output <comparison.json>
```

The analyzer requires consecutive factor-two grids with the same time step,
endpoint, initial profile, mechanism and native energy convention. It verifies
complete snapshots and analytic fluid fractions, restricts conserved densities
with fluid-volume weights, and checks the integral before and after restriction.
Covered coarse cells never contribute. Stored temperature is volume averaged
for comparison; it is not reconstructed from the restricted conserved state.
Initial, endpoint and evolution-increment RMS/max differences are separate.
With three grids the output includes observed orders; two grids report only
differences. A fixed nonzero time step and three grid sizes alone do not
establish the asymptotic spatial order or physical validation.

To separate time-step effects from spatial differences, generate two or three
study captures on the same mesh with factor-two root time steps and identical
endpoint time, then use:

```sh
python3 tools/compare_tgreact_eb_meshes.py --temporal \
  --runs <dt-run> <half-dt-run> <quarter-dt-run> --output <temporal.json>
```

This mode compares corresponding fluid cells without spatial restriction.
Initial fields must match exactly; profile, chemistry, EOS and transport
controls must also match. Each trajectory initializes its own reaction history
for its step size. Three captures yield observed temporal orders, with null
orders for roundoff-level differences. These diagnostic orders do not qualify
the original reference or establish asymptotic convergence by themselves.

Mesh-study JSON also accepts `initial_velocity_cm_s: 0` for a separate
initially stationary control. It overrides native `prob.u0`; PeleF imports
that same zero-momentum initial state, including its consistent total energy.
Reaction, transport, EB geometry and AMR remain enabled. The checker requires
exactly zero initial speed in this mode, while the default retains the frozen
minimum-speed requirement. A separate weak-flow control accepts 1 cm/s,
one thousandth of the original 1000 cm/s amplitude, and scales only the
initial minimum-speed admission by that ratio. Field and conservation limits
remain unchanged. These three amplitudes are the only accepted values.
Reports record the control and the mesh/time analyzer refuses to
mix amplitudes. This intervention changes the initial flow, not the hydro
operator; pressure-driven motion can still develop and conclusions about
limiters or interfaces require additional evidence.
