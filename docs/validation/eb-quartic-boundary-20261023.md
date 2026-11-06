# Runtime attribution at the remaining shear boundary

This is isolated serial 16-cell/50-ns candidate evidence, not C2 completion.
Artifacts: `build/eb-quartic-boundary-20261023/`. Build commands, source identity,
frozen candidate linkage and binary SHA are in `build-record.json`; execution
and four-table identity checks are in `grid16-record.json`.

The diagnostic run completed successfully. All four ordinary tables match the
saved candidate byte-for-byte; historical stage and sync rows are unchanged.
The runtime parent (16,1,8) has all three stencil directions available but fails
the smoothness screen on all four recorded fills. Its initial y momentum second
differences are [0.111940707, -0.111940707, -0.318780163], versus roundoff bound
1.11334e-14: this is a resolved curvature-sign change, not the earlier near-zero
transverse roundoff problem. `admitted=true, theta=1` here describes the fallback
family; it does not prove the rejected high-order family passes EOS checks.

At fine halo (0,1,15), the initial runtime velocity is 1.006083715962 m/s,
matching the analytic MC linear interpolation within 1e-12. The exact fine-cell
average is 0.978597633264 m/s, leaving +0.027486082698 m/s. `audit_boundary.py`
checks these facts from the measured halo and prepared donor lines.

Applying the pinned quartic coefficients to that saved initial y line gives
0.979470818621 m/s, with error +0.000873185357 m/s, or 3.17683% of the linear
error. `replay_initial.py` records this counterfactual interpolation only.
It does not establish weighted-family EOS validity, later-stage admission,
trajectory improvement, discontinuity robustness or MPI/restart behavior.

Combined with the stage attribution in the mesh record, this identifies an
actual fallback-induced initial velocity mismatch at the remaining maximum-error
boundary. A bounded causal intervention is still required to quantify how much
of the later pressure/momentum error this mismatch explains. No production
predicate or numerical method was changed. Python syntax, project contract and
`git diff --check` pass.
