# Attribution of the reacting halo-policy difference, 2026-10-22

The large change from the pinned reacting trajectory is already present in
coarse-to-fine interpolation before chemistry or time integration. The
separate-momentum option changes both scalar and momentum slopes by design.
The recorded kernel agrees with an independent evaluation of that design;
this audit found no stencil-indexing or limiter-formula implementation error.

## Method and identity

The exact state records from the completed 8-cubed reacting comparison are
reused. The first attempt with its V1 header was refused by the existing
initial-audit admission check. A separate V4 diagnostic input enables the
converged EOS and tight reaction tolerances without changing any state,
geometry, spacing or step row. No reaction or evolution is performed here.
The changed EOS admission applies to pressure observations; the limiter
operates directly on the unchanged recorded vectors and is independent of
that EOS setting. The V1 refusal log is retained, not counted as a pass.

Artifacts: `build/eb-momentum-attribution-20261022/`. The frozen probe and
`record.json` identify both initial-operator executions, exact arguments,
input/binary/source hashes and durations. `initial_limiter.csv` files record
center, six neighbors, offsets, and actual default/separate kernel outputs
for every component on each eligible fluid coarse/fine face stencil.
`tools/analyze_eb_halo_limiter.py` independently evaluates the scalar and
momentum limiter factors and compares its reconstruction to both outputs.
It accepts the LiDryer layout explicitly: 14 gas components, stored
temperature, then held internal energy.

## Results

- 128 fluid face stencils, 384 directional limiter evaluations.
- 48 directions change the shared scalar factor when momentum is excluded.
- In 24 directions, the default factor is zero while the scalar-only factor
  is nonzero. Resolved scalar slopes are therefore restored by the option.
- Independent kernel reconstruction matches within 1.12e-16 scaled error.
- Both policies reproduce actual production halos exactly in the audit;
  symmetric-child mean residual is at most 2.49e-16.

For example, the fine ghost at (5,0,7) uses common factor zero in x because
of the x-momentum limiter, whereas the scalar-only factor is one. The
centered x-momentum slope is 0.21819448 kg/m2/s, a resolved variation rather
than a roundoff-sized perturbation. At that ghost the packed temperature
changes from 1980.865828 to 1989.031847 K (8.166019 K), density changes by
6.63892e-4 kg/m3 and held-energy density by 357.3195 J/m3. This explains why
separate momentum is a broad numerical-policy change in a Taylor–Green flow,
not merely a rest-state roundoff repair.

| Initial-operator observable | Default | Separate momentum |
| --- | ---: | ---: |
| Maximum adjacent ghost pressure offset (Pa) | 0.0366603 | 0.275222 |
| Maximum area-weighted coarse/fine mass-flux mismatch (kg/m2/s) | 0.181588 | 0.0274046 |
| Maximum normal-momentum flux mismatch (Pa) | 116.546 | 35.5533 |
| Maximum energy-flux mismatch (W/m2) | 596276 | 89673.2 |

Open coarse/fine areas match exactly. These are instantaneous intensive
flux mismatches, not accumulated conservation errors. This variable-velocity
flow has physical flux; zero mass flux is not an acceptance condition.
The option reduces the flux mismatches while increasing the initial pressure
offset. Neither that improvement nor unchanged global budgets establishes
better continuum accuracy. Conversely, a larger difference from the pinned
trajectory does not by itself prove worse continuum accuracy.

## Decision and verification boundary

Do not restore shared momentum/scalar clipping just to force the new policy
back toward the pinned trajectory: that would restore the known constant-
velocity/rest-state defect. No production limiter or default is changed here.
The pinned default path still owns its original field-comparison gate. The
opt-in policy needs independent variable-velocity, all-field convergence
against a resolved or analytic reference; its constant-velocity pressure
mesh result is insufficient for that purpose.

The added probe output is diagnostic only. Build, both initial-operator
runs, independent stencil replay, mechanism provenance, Python syntax,
project contract and whitespace checks pass. The existing audit outputs
are compared byte for byte against the preceding frozen probe; the result
is in `parity.json`. This is not a new reacting-trajectory or MPI qualification.
