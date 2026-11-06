# Optional coarse-halo momentum policy: application gate, 2026-10-21

The selected H2/O2/N2 gas EB-tree application passes common and recursive
held-MOL lifecycle checks with `separate_momentum=true`. This is an opt-in
numerical policy; it does not close original-PeleC parity, spatial/temporal
convergence or physical validation. Default remains false. Spray, enabled
LES and prepared redistribution are excluded.

The three independently rerun CTests passed with zero failures or skips:
default converged-EOS CVODE restart, candidate common restart and candidate
recursive restart. Test times were 46.52, 54.39 and 92.07 s (193.02 s total).
The existing isolated Debug build uses GNU Fortran 13.3 and OpenMPI 4.1.6,
with the configured SUNDIALS backend and pinned h2o2_full bundle. Build,
provenance, Python syntax, project contract and whitespace checks passed.
The MPI runtime dependency audit passed (ABI 40). Compiler warnings for
unused procedures and the existing executable-stack requirement remain.

The reacting test includes molecular transport, topology collapse/recreation,
load redistribution, serial/MPI fields and snapshots, reaction-off control,
and two-rank checkpoint continuation on three ranks. Both directions of
checkpoint policy mismatch and mixed-rank input mismatch are refused before
result output. Checkpoint bytes remain unchanged after rejection.

| Clock | Serial maximum budget residual | Serial hydrogen reacted fraction |
| --- | ---: | ---: |
| Common | 7.3039e-16 | 2.2952255e-4 |
| Recursive | 8.2268e-16 | 2.2952263e-4 |

The endpoint is 300 ns. These are bounded application inventories, not an
ignition delay or a reference-solver comparison. Serial chemical-element
residuals are below 6e-16. Fixed-mechanism coverage also passes: the original default CTest and fresh
common/recursive candidate lifecycle runs. Candidate inventory residuals are
`1.097e-15` and `1.132e-15`. Both restart mismatch directions, serial/MPI
unsupported configurations and late invalid-restart output preservation pass.
The first candidate runs reached the legacy non-held smoke with the candidate
flag still enabled and were correctly refused. The fixture now disables that
flag with held energy; no numerical source or tolerance was changed. Original
failure logs remain in `build/eb-momentum-fixed-application-20261021/`; corrected
commands, hashes and fresh outputs are in its `corrected/` directory.

Exact source hashes, commands, build/test logs and return codes are in
`build/eb-momentum-application-20261021/record.json`. Per-run input/executable
hashes and checks are linked by `verified-reports.json` in that directory.
No commit identity alone describes this dirty checkout; use these hashes.
Kernel, operator, retry and complete-window evidence is linked from the
[policy design](../design_decisions/eb-halo-momentum-limiter.md).

## Default-path regression closure

Seven additional existing CTests passed after the optional interface changes:
serial full-tree integration, serial/MPI MOL agreement, rotating MOL,
recursive and local interval histories, and common/recursive StateRedist.
There were zero failures or skips; total test time was293.49 s. Source hashes,
commands and logs are in `build/eb-momentum-default-regression-20261021/`.
This is focused changed-path regression, not the full configuration matrix
or an installed-binary audit. No further solver job is pending for this
implementation increment.
