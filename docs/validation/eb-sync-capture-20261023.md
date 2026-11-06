# Recursive held-energy closure capture — 2026-10-23

The serial recursive SDC path now accepts an optional read-only
`observe_mol_sync` callback around the ordinary held-energy closure. Phase 0
precedes the closure; phase 1 is emitted only if it succeeds. The callback
receives gas, held internal energy, integrated gas/scalar face registers,
level, interval start and interval width. This measures the combined closure
(reflux, restriction and energy reset), not each internal suboperation.
It covers the non-redistributed ordinary held path; other closure/replay paths
are not claimed. As with stage capture, observations are attempts, not
accepted-step receipts, and times are relative to the recursive invocation.

The shear stage probe records five leaf locations: the two fine maxima probes,
two coarse cells across their periodic x boundary, and one coarse cell above
the fine x/z corner. The output stores six integrated face registers per cell
and pre/post momentum, held energy and EOS pressure. These flux units differ
from the instantaneous spatial-stage CSV; no division by dt is applied.

The bounded 16-grid, 50 ns cut-shear comparison will test unchanged ordinary
outputs and unchanged spatial-stage traces, then compare the closure momentum
increment with the divergence of changed face registers. This separates an
additional closure contribution from the already observed initial transverse
acceleration. No numerical policy is changed.

Artifacts: `build/eb-sync-capture-20261023/`. Source backups precede this change.
The build and cut-shear probe completed successfully (54.73 seconds).
All 16 ordinary tables and all four spatial-stage traces are unchanged.
The initial independent closure audit **failed**: it assumed that every leaf
momentum change was solely the divergence of changed face registers. At the
root closure, the fine cut probe has no changed face registers but its
x-momentum changes by 1.325e-9 (default) or -1.218e-7 (paired). This is not
roundoff at the existing gate, and the gate has not been relaxed.

Inspection of `close_held_prepared_eb_tree_3d` shows that even this
`redistribution=none` path redistributes the reflux mismatch when extended
geometry is supplied and a finer level exists (`flux_rereflux`). It forms
joint gas/held-energy increments, constructs re-reflux registers, then routes
an additional cell correction. The face registers alone therefore omit that
part of the update. The initial failed audit, record and CSVs are retained
unchanged; they do not establish a numerical conservation defect.

At the two fine probes, this first 50 ns root closure changes neither y nor
z momentum. The coarse cell above the fine corner receives a z-momentum
increment -7.74307e-8 matching its changed face registers. These observations
support separating the already present fine transverse acceleration from
closure effects, but do not exclude later-time feedback through coarse halos.

The follow-up adds optional raw and actually applied gas increments to the
closure diagnostic. It will check (a) changed face registers against the raw
increment, (b) applied increment against the observed leaf update, and (c)
the complete two-stage fine momentum trajectory including all measured
closure changes. Artifacts are `build/eb-sync-increment-20261023/`.
This corrects the diagnostic accounting rather than altering the solver.

The first build of the callback rejected a repeated `present` test of a
host-associated optional procedure in the recursive worker; evaluating presence
once in the outer routine resolved it. `build-first-failed.log` retains that
compiler failure. Mechanism provenance, project contract, Python syntax and
whitespace checks passed. No MPI or fixed-PeleC qualification is claimed.
