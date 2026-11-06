# Fixed-mesh reacting EB-AMR time refinement, 2026-10-21

All four new trajectories pass the unchanged field and budget gates. The
8- and 16-cubed studies compare 50, 25 and 12.5 ns root steps at 100 ns,
reusing the completed 50 ns smooth-profile captures. Both solvers retain the
native activation convention and held spatial EOS; PeleF uses its explicit
converged caloric inverse. Chemistry tolerances remain relative `1e-12` and
absolute density `1e-18 kg/m3`. These are diagnostic comparisons, not
qualification of unmodified PeleC or the canonical mechanism.

## Observed temporal orders

The native held-EOS orders agree with these PeleF values to the shown digits.

| Field | 8 cubed | 16 cubed |
| --- | ---: | ---: |
| Density | 1.061 | 1.111 |
| x-momentum | 1.098 | 1.095 |
| z-momentum | 1.108 | 0.929 |
| Total energy | 1.071 | 1.113 |
| Temperature | 1.079 | 1.028 |
| H2O density | 1.339 | 1.566 |

These are volume-weighted endpoint self-differences on fixed meshes. Initial
states match exactly within each solver. There is no evidence here of
second-order temporal convergence or of a Fortran-only temporal defect.

## Time effects versus spatial differences

Conservative 16-to-8 restriction at each common time step gives these
PeleF RMS differences of the evolution increments, `state(100 ns)-state(0)`:

| Field | 50 ns | 25 ns | 12.5 ns |
| --- | ---: | ---: | ---: |
| Density, kg/m3 | 6.57232e-6 | 6.54118e-6 | 6.52767e-6 |
| x-momentum, kg/m2/s | 3.25504e-3 | 3.23742e-3 | 3.23083e-3 |
| z-momentum, kg/m2/s | 2.62010e-4 | 2.71859e-4 | 2.76078e-4 |
| Total energy, J/m3 | 22.1046 | 21.9982 | 21.9520 |
| Temperature, K | 0.0400973 | 0.0398940 | 0.0398052 |

Density, total energy and temperature spatial differences change by less
than 1% when the step is quartered. The z-momentum difference instead grows
by about 5.4%. Thus time-step effects alone do not explain the observed
spatial differences on this pair. This does not establish a spatial order
on finer meshes or isolate limiter, interface and boundary contributions.

The next separate control sets the initial TG amplitude to zero while retaining
hydro, reaction, transport and EB-AMR. It tests dependence on the initial flow;
it is not equivalent to switching off hydro or proving limiter causality.

## Implementation and reproducibility

`tools/compare_tgreact_eb_meshes.py --temporal` now compares completed fixed-mesh
captures with factor-two time steps and a common endpoint. It rejects mixed
physical/numerical controls and differing initial states, validates complete
snapshots, and retains null orders for roundoff-level differences. Existing
spatial outputs are unchanged: all solver results from the saved three-grid
smooth-profile study were reproduced exactly.

The three Python tests pass, including a manufactured second-order temporal
sequence and mixed-control refusals. Python syntax, H2/O2 and LiDryer provenance,
project contract and `git diff --check` pass. These are checker and diagnostic
checks, not a fresh MPI, installed-release or full configuration matrix.

Artifacts: `build/eb-smooth-temporal-20261021/`. `record.json` contains exact
commands, durations, report paths and hashes; `source.diff` preserves the
tracked working-tree delta at launch. HEAD is
`0a0e00ebe9138145fbb0e0e745665d0d768df009`, with existing uncommitted work.
The reused Debug probe and its link inputs are identified by
`build/eb-eos-convergence-20261021/probe-build.json`; they were hash-checked
before execution. This study did not rebuild the solver. Checker copies are
preserved as `*_at_run.py`. `analysis.json` joins the saved time and spatial
results. Wrapper exit is zero; elapsed wall time is 1478.71 s (24.65 min).

The original-reference field gate, asymptotic convergence, broad release
qualification and physical validation remain open.
