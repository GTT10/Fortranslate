# Initially stationary EB-AMR control, 2026-10-21

The 8- and 16-cubed stationary controls pass the unchanged trajectory gates.
The 32-cubed solvers both finish, but their comparison fails 19 field limits.
No three-grid order is reported: the mesh analyzer requires passing captures.
This is a newly exposed diagnostic discrepancy, not a timeout or successful
fixed-reference qualification.

The control sets native TG `prob.u0=0` and imports the resulting consistent
gas/energy state into PeleF. Every initial momentum is zero. Chemistry,
hydro, molecular transport, EB and recursive AMR remain enabled. Other
settings match the smooth129, native-activation/held-EOS comparison: 50 ns
steps to 100 ns, relative chemistry tolerance `1e-12`, absolute density
tolerance `1e-18 kg/m3`, and explicit converged PeleF caloric inversion.
Zero initial speed is checked exactly; the original nonzero-flow case is
unchanged.

| Grid | Result | Maximum endpoint temperature difference |
| --- | --- | ---: |
| 8 cubed | Pass | 1.0961e-5 K |
| 16 cubed | Pass | 7.8479e-6 K |
| 32 cubed | Fail, 19 field limits over two times | 9.3046e-3 K |

On the 32 grid, maximum density differences are `1.5778e-7 kg/m3` at 50 ns
and `4.4254e-7 kg/m3` at 100 ns. The initial temperature difference is zero;
at 50 ns it reaches `7.8679e-3 K`. Density and temperature maxima occur on
the fine patch's upper boundary. Energy reaches `1.6396 J/m3` difference at
100 ns. Species/element inventory and energy-budget gates do not fail;
conservation does not substitute for field agreement.

## Saved-stencil investigation

`build/eb-stationary-halo-audit-20261021/audit_scales.py` evaluates the shared
coarse-halo slope limiter on saved coarse states, recording input hashes and
the components controlling each coefficient. It follows
`interpolate_eb_linear_halo_3d` and the pinned AMReX
`mf_cell_cons_lin_interp_llslope`: the minimum component ratio multiplies
every component's slope in a direction.

Initial coefficients agree on all three grids. At 50 ns, transverse momentum
slopes of order `1e-33` can force the common coefficient to zero in one
implementation while the other retains one. For example, the 32-grid coarse
cell `(9,4,18)` in the x direction has native coefficient zero controlled
by y-momentum and PeleF coefficient one. The corresponding difference in
density slope is about `7.8052e-4 kg/m3`. The SI unit-conversion control does
not explain that example. Similar sensitivities occur on the passing smaller
grids, so their presence alone does not explain why the 32-grid gate fails.

This is evidence of a discontinuous response near exact rest in the shared
limiter. These are saved endpoint stencils, not captured coarse-time donors
at the actual fine-stage evaluations. They identify a causal candidate, not
proof that this branch causes every observed field failure. No tolerance,
production limiter or reference hash has been changed.

A separate 1 cm/s TG amplitude control, one thousandth of the usual value,
is implemented to test dependence on exact rest. Its initial minimum-speed
admission scales with the amplitude; field and conservation limits do not.
It runs the previously failing 32 grid first and starts 8/16 only after a
passing result. It remains an explicitly identified diagnostic.

The weak-flow 32-grid control has now completed. At 50 ns every field gate
passes; at 100 ns three limits remain: scaled density `1.6616e-8`, total
energy `7.1904e-8`, and held energy `7.1903e-8`, against `1e-8` limits.
Endpoint temperature difference falls to `9.0715e-5 K`, below its `1e-3 K`
limit. The driver correctly stops before 8/16 because 32 did not pass.
This establishes strong dependence on the initial-flow amplitude, but does
not establish convergence or prove the shared limiter is the sole cause.

The next diagnostic captures live fine ghost values and selected coarse
limiter stencils during the first 50 ns of the stationary32 failure. It uses
private copied sources/objects and leaves production numerical operations
unchanged. Native and Fortran observer endpoints must be byte-identical to
their respective saved endpoints before the trace can be interpreted.
Artifacts are in `build/eb-live-halo-trace-20261021/`. The first native
observer compile failed on a diagnostic-only component-count symbol; using
`S.nComp()` fixed it. The resumed native and Fortran runs finish, and their
50 ns endpoint exports match the original respective exports byte-for-byte.
Initial fine ghost inputs also agree exactly. However, the native trace cap
counted internal tiles, covering only the first three fine spatial evaluations.
A corrected native-only observer records full ghosts once per spatial
evaluation, including time/dt/reflux metadata. That capture completes in
`native-full/`, with another byte-identical native endpoint; all twelve
Fortran halo captures are reused. Six native spatial evaluations correspond
to the two fine substeps, each with old/predicted/corrected evaluations.
Fortran records two prepared halos per spatial evaluation. All 5584 recorded
fluid ghost cells have matching index sets. Initially conserved inputs agree
to roundoff; the interpolated native temperature and EOS-recovered Fortran
temperature are distinct quantities and are not treated as a separate field
failure in this ghost audit.

### Live local intervention

At the first fine predictor evaluation, native call 2 and Fortran calls 3/4
use the half-parent-time donors. At fine ghost `(17,10,33)`, coarse donor
`(9,5,17)`, Fortran's x-direction factor is one but native interpolation is
constant. Conservative eight-child averaging of the native fine ghost data
recovers the coarse y-momentum stencil: `-7.0631e-26`, zero, and
`-7.0513e-31 kg/m2/s`. Fortran records exact zeros in those three entries.
The native tiny local maximum forces the shared x factor to zero; the exact
zero Fortran stencil imposes no constraint.

Replacing only those three y-momentum values in the saved live Fortran
stencil changes its factor from one to zero. Density ghost error drops from
`1.95130e-4` to `6.8001e-16 kg/m3`; energy error drops from `105.535` to
`4.1910e-9 J/m3`. The reconstructed original Fortran interpolated density
also matches its captured value before intervention. This establishes the
local amplification mechanism, not yet the cause of the full endpoint error.
The same shared-factor rule exists in both implementations; injecting native
roundoff into production or widening field tolerances would not be a fix.

`local_intervention.py` and `local-intervention.json` preserve this calculation.
A matched diagnostic ablation is now running in
`build/eb-limiter-ablation-20261021/`: copied native and Fortran interpolators
exclude momentum from selecting the shared factor while retaining common
scaling of all component slopes. It tests the failing first 50 ns with the
original field and conservation gates. This is not an adopted interpolation
policy or canonical-reference qualification. Its first setup attempt used
the wrong archive path and stopped before compiling; the corrected attempt
uses the selected regular-3D runtime that owns the transfer module.

The matched momentum-factor ablation has completed successfully at 50 ns.
Every unchanged field and budget gate passes. Maximum differences are
`4.2674e-14 kg/m3` in density, `1.6365e-7 J/m3` in total energy, and
`1.1401e-5 K` in temperature; the relative energy-budget residual is
`8.1987e-17`. This extends the local intervention to the failed first-step
trajectory and identifies shared momentum limiting as the amplification
mechanism in this diagnostic. It does not qualify the original solver or
prove an unrestricted momentum slope is safe.

A separately bounded momentum-limiter kernel is now implemented as an opt-in
candidate; application defaults remain unchanged. Its design and promotion
gates are in [the policy note](../design_decisions/eb-halo-momentum-limiter.md).

## Evidence and checks

Stationary artifacts: `build/eb-zero-velocity-mesh-20261021/`. `record.json`
contains all commands, source hashes, per-grid reports and costs. The wrapper
exit is 1 after 2575.01 s (42.92 min); native and Fortran solver exits are
zero. The failing `checks.json` retains all 19 violations and all three
snapshots. The previously qualified Debug probe is reused with binary/link
hash checks; this is not a fresh compiler/MPI/install audit.

Saved-stencil evidence: `build/eb-stationary-halo-audit-20261021/scales.json`.
Weak-flow job records: `build/eb-weak-velocity-mesh-20261021/`; wrapper exit 1
after 2188.84 s (36.48 min), with both solver exits zero and three failed
field limits. Four Python tests pass, including bounded-velocity
admission/refusal and mixed-amplitude comparison refusal. Python syntax,
LiDryer provenance, project contract and `git diff --check` pass.

Canonical-reference agreement, asymptotic convergence and physical validation
remain open. In particular, the failed stationary fine grid is retained as
an unresolved case rather than removed from the acceptance scope.
