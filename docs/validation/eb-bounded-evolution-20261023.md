# Bounded tensor candidate: completed 16-cell evolution

The isolated serial paired-policy 16-cell/50-ns candidate completes in
54.483 seconds with normalized budget residual 5.92935e-14. All six tables
are complete and finite, and initial analytic rows match the saved baseline.
Stage directional face-rate closure, 27 synchronization pairs and six fine
momentum trajectories pass, including wall and measured closure increments.

Against the saved vector-roundoff quartic candidate, final y-momentum
L1/L2/Linf ratios are 0.0318824/0.0302425/0.0317692; z ratios are
0.0343011/0.0320659/0.0317673; pressure ratios are
0.0342674/0.0320577/0.0317673. Final maxima are y=2.70956e-7,
z=2.77099e-7 (momentum density units) and pressure=0.00681394 Pa.
The four recorded interface probes have first-predictor pressure error
0.003546807 Pa, versus 0.111649534 Pa in that baseline.

No norm increases among the 17 components and three norms versus either
that baseline or the separate inflection-hierarchy candidate. This removes
the latter candidate's measured y-L2 regression at this condition. It does
not establish fixed-PeleC parity: both baselines are isolated diagnostic
variants. Multigrid convergence, longer-time behavior and production/MPI
integration remain unverified; C2-C5 remain open.

The frozen tree changes only the call to the isolated family routine,
enabling bounded_correction for all eligible cached donor families. The
existing geometric scope, volume guard, EOS admission and fallback are
retained. The tensor averages all six directional orders; measured runtime
is about 9% greater than the saved 49.87-second inflection run, but these
single wall-clock observations are not a controlled performance benchmark.

Source/binary identity and compiler/link commands are in
`build/eb-bounded-evolution-20261023/build-record.json`; the completed receipt
is `grid16-record.json`. The input SHA is
`1d8f68616256d990eaacaa6bfa2ea1c9f61115ac5b8a95d3930507b0cc9ddc2e`.
`grid16-audit.json`, `grid16-trajectory-audit.json`, `comparison.json` and
`all-comparison.json` retain the quantitative evidence.

Analysis commands (after the completed run):

```sh
python3 build/eb-bounded-evolution-20261023/audit.py --grid 16
python3 build/eb-bounded-evolution-20261023/summarize_forcing.py --grid 16
python3 build/eb-bounded-evolution-20261023/compare.py
python3 build/eb-bounded-evolution-20261023/compare_all.py
```

## Completed 32-cell refinement discriminator

The same executable completed the 32-cell/50-ns run in 366.608 seconds,
return code zero. All six output tables are complete and finite; maximum
normalized budget residual is 7.58551e-14. Stage face closure, all 27 sync
pairs and all six fine momentum trajectories pass. Compared with the saved
32-cell vector-roundoff quartic candidate, none of the 51 field norms
increases. The y/z momentum and pressure Linf ratios are respectively
0.00774890, 0.00801232 and 0.00801756.

| Quantity | 16-cell Linf | 32-cell Linf | Observed log2 ratio |
|---|---:|---:|---:|
| y momentum density | 2.70955959e-7 | 3.30171544e-8 | 3.03677 |
| z momentum density | 2.77099257e-7 | 3.49158856e-8 | 2.98845 |
| pressure (Pa) | 0.00681393924 | 0.000422986253 | 4.00981 |

No norm increases on refinement. The comparison checks the same binary,
physical domain, plane, EOS/stepping controls, zero initial analytic errors,
50-ns endpoints and successful stage/trajectory audits. Initial string-based
control comparison rejected the equivalent `5e-8` and long decimal spelling;
it was corrected to exact parsed numerical equality, with no tolerance change.
The initial failed analysis is not a solver failure or a loosened accuracy gate.

These are observed ratios on two grids, not demonstrated asymptotic orders.
Short-time forcing and localization can produce high integrated-norm ratios;
no global fourth/fifth-order scheme claim follows. The original coarse-grid
scalar limitation remains, and neither comparison uses fixed PeleC as its
baseline. Production/MPI/restart integration and longer-time behavior remain
open. Next use a matched regular-geometry control and saved longer-time
conditions before promoting the candidate; do not launch a larger grid merely
to obtain another favorable slope.

Additional artifacts: `grid32-record.json`, `grid32-audit.json`,
`grid32-trajectory-audit.json`, `comparison32.json`, `refinement.json`.
The refinement checker preserves all 51 norms, both-grid maximum locations,
and table/binary identities. Python syntax, project contract and whitespace
checks pass. No job remains active.
