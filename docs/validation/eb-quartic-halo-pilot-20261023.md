# Frozen quartic halo pilot — 2026-10-23

## Decision

The first initial-operator pilot exceeded its 180-second candidate limit.
Do not count it as a completed four-policy comparison or launch evolution
from this result. The baseline completed in 16.6995 seconds. No production
default, checkpoint policy or reference tolerance changed.

## Scope and source identity

Artifacts are in `build/eb-quartic-halo-pilot-20261023/`. `build-record.json`
records compile/link commands and hashes of the frozen tree-module copy,
helper sources and executable. `initial-record.json` records both executable
hashes, input hash, commands, process identities and the timeout. The driver
terminated only its own candidate subprocess. The job handle subsequently
returned exit 0 because the driver records failure in JSON; **the comparison
itself failed**, as `status=failed` and `timed_out=true` show.

The candidate replaces only the copied tree module, for a two-level root-donor
case with periodic x/y and slip-wall z boundaries. It gathers a five-cubed
stencil using existing ancestor-aware state/held-energy sampling, applies the
connected-geometry and weighted-family checks, and adopts only smooth,
positive-theta families. It is not a general MPI halo implementation or the
pinned PeleC interpolation policy. This initial audit uses the existing four
policies without an evolved connected-directions comparison.

## Salvaged evidence

`partial-audit.json` compares available tables with the fresh baseline. Default,
momentum and velocity policies each emitted pressure (22,736 rows), flux
(26,112 rows) and rates (58,240 rows), with matching row counts. No paired-policy
result or successful final receipt is available from the interrupted run.

- All three default tables are byte-identical to the baseline.
- At fine ghosts (0,8,11) and (0,9,16), momentum/velocity policies change u
  from 9.7449535840 to 9.9356948241 m/s. Relative to the corresponding analytic
  fine average, the local gap decreases from about 0.190915 to 0.000173687 m/s.
- Across the recorded domain, maximum absolute density, streamwise momentum
  and energy initial rates remain 584.8697, 352795.3548 and 381320788.6034 in
  their respective SI rate units. Transverse initial rates remain zero.
  A local ghost improvement therefore does not establish global operator or
  evolved maximum-error improvement.

`partial-localization.json` and its reproducer `localize_partial.py` separate
coarse and fine rates. All recorded coarse rates are zero. Fine density has
its maximum at (1,2,12), while streamwise momentum and energy maxima are at
(1,2,11); these are different from the two selected transverse-error probes.
At (1,8,11), the density rate changes from -117.770434 to -0.107159949,
streamwise momentum from -71868.393857 to -65.393601, and energy from
-77313488.0566 to -70354.0477. The (1,9,16) rates show the same reduction.
Thus the local initial forcing is reduced, but the global forcing maximum
persists elsewhere along the interface. This must be followed through pressure
generation and subsequent transverse rates before claiming an evolved gain.

The next diagnostic copy caches a family only within one halo-fill call and
only for the same root donor index, child-volume weights and temperature guess.
It preserves rejection as well as adoption. This avoids recomputing all eight
children for each sibling ghost. Its results must match all nine available
uncached tables before its runtime or new paired output can be accepted.
The cache variant is in `build/eb-quartic-halo-cache-20261023/`; it remains
outside production source and CMake.

## Completed cache pilot

The cache copy built successfully and completed all four initial policies in
**38.5896 seconds**, within the same 180-second limit. All nine available
uncached CSV tables match exactly; default also matches the baseline. For each
of momentum, velocity and paired pressure, the candidate reports 960 adopted
ghosts and 6,888 cache hits. Cache hits include rejected families, so these
counts have different meanings. The fresh paired results also retain the
global initial maxima above. `initial-record.json`, `initial-audit.json`,
`candidate.log` and the frozen build record preserve the evidence in the
cache artifact directory.

A bounded 16-cell/50-ns evolution with stage and closure traces is now launched
in `build/eb-quartic-evolution-20261023/`, using this identical candidate binary
and the saved connected-direction baseline input. The run limit is 240 seconds;
the driver records failures rather than admitting incomplete comparisons.
The audit requires unchanged default outputs and independently closes face
rates, raw/applied synchronization increments and fine momentum trajectories.
Candidate norms and maximum locations must be compared separately. Evolution
and convergence conclusions remain pending.
