# Reacting paired-pressure comparison against the fixed reference, 2026-10-22

The paired pressure policy was compared against an unchanged captured pinned
PeleC TGReact EB-AMR trajectory. This is a reacting and molecular-transport
comparison, separate from the completed hydro-only pressure mesh diagnostic.
It does not replace that reference or change any acceptance threshold.

## Case and provenance

Two 8-cubed patches, ratio two, a 1 cm domain, planar EB at z=3.3 mm,
Taylor–Green velocity, nonuniform temperature and H2 composition, three
100 ns recursive MOL2 root steps, LiDryer chemistry, and Simple transport.
Initial reaction history is computed independently for each run. The frozen
33-point profile and original EOS/reaction controls are unchanged.

Capture: `build/tgreact-eb-amr-reference-20261019/runs/matched_eb_amr_62p3fjsx/reference`.
The checker verifies every file in the capture manifest, pinned source and
input contracts, and transport-fit provenance before reuse. No new native
reference simulation was needed. Current artifacts are under
`build/eb-pressure-reacting-20261022/`; `record.json` records exact commands,
source and binary hashes, reference identity, durations and report paths.

The probe now accepts `--paired-reacting`; the checker exposes
`--paired-pressure`. These enable separate momentum, consistent velocity and
paired pressure in the ordinary reacting/transport driver. A successful
intervention comparison would still be marked `qualification=false` because
it changes the pinned spatial algorithm. The ordinary two-argument probe
keeps the previous policy. A run with the previous frozen executable matches
all eight default trajectory/budget CSVs byte for byte.

## Completed fixed-reference comparison

Both simulations complete all three steps. Both fail the unchanged field
comparison gates; these are numerical comparison failures, not execution
failures. `comparison.json` retains every field and failed gate. At 300 ns:

| Observable and metric | Default policy | Paired policy | Fixed limit |
| --- | ---: | ---: | ---: |
| Temperature, max absolute K | 1.09800e-4 | 2.38882e-1 | 1e-3 |
| Density, max difference/reference max | 5.74765e-8 | 3.99888e-4 | 1e-8 |
| Total-energy density, same relative metric | 8.48018e-8 | 5.60069e-4 | 1e-8 |
| x momentum, max absolute kg/m2/s | 4.76496e-6 | 7.05729e-2 | 1e-6 |
| H2 mass fraction, max absolute | 5.06984e-10 | 1.34430e-6 | 1e-7 |

Correction-accounted relative energy residuals are at most 8.09e-17 (default)
and 9.47e-17 (paired). Relative H2-consumption errors are 2.242e-5 and
2.251e-5, both within 1e-3. Conservation and chemical-consumption agreement
do not repair the spatial field failure. The policy improves the separate
constant-velocity pressure diagnostic but worsens agreement with this
variable-velocity reacting reference. This alone does not determine which
trajectory is closer to a continuum solution.

## Attribution

All four policy configurations completed the same three-step trajectory;
all four fail at least one unchanged field gate. The two additional runs
use `--momentum-only` and `--velocity-only`, forwarded to explicit probe
options. `ablation/record.json` records the separate rebuilt executable and
source identities; `ablation-comparison.json` records all four results.

| Enabled policy | Temperature error (K) | Density relative error | x momentum absolute error (kg/m2/s) |
| --- | ---: | ---: | ---: |
| Default | 1.09800e-4 | 5.74765e-8 | 4.76496e-6 |
| Separate momentum | 2.40975e-1 | 4.02108e-4 | 7.08118e-2 |
| Plus consistent velocity | 2.40260e-1 | 4.01762e-4 | 7.10491e-2 |
| Plus paired pressure | 2.38882e-1 | 3.99888e-4 | 7.05729e-2 |

The dominant departure from the pinned trajectory is already introduced by
separate momentum. Consistent velocity and paired pressure do not explain
the jump from the default trajectory. Paired pressure reduces temperature
error by about 0.57% relative to consistent velocity alone; this is not a
passing reference result. The shared material/temperature limiter and its
coupling to momentum therefore need a variable-velocity accuracy assessment,
not a pressure-only repair chosen to hide the discrepancy.

Across all four runs the maximum correction-accounted relative energy
residual is 1.22e-16. Every run completes but comparison acceptance fails;
none is counted as a passed fixed-reference trajectory. No production
numerical algorithm, default, tolerance or golden hash was changed in
response. There are no remaining jobs in this comparison batch.

## Verification

The probe builds successfully with GNU Fortran 13.3.0 in the existing Debug
OpenMPI configuration. The new command-line policies exercise the actual
reacting core; no private time integrator is substituted. LiDryer and full
H2/O2 mechanism provenance, Python syntax, project contract and
`git diff --check` pass. The default eight-CSV byte comparison was completed
for the first paired-option build against the preceding frozen executable;
the subsequent ablation build adds the two explicit prerequisite selectors.
