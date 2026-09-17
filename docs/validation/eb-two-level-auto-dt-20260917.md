# Two-level EB application automatic interval — 2026-09-17

Implementation: `112eacb3f72bbcccb25c8bbe81911cd9690915e7`, based on
`c401950` on local `main`.

The fixed-mechanism serial application now accepts `automatic_timestep`
(default false), `hydro_cfl` (default 0.5, allowed `(0,1]`) and `transport_cfl`
(default 0.35, allowed `(0,0.5]`). When enabled, each step privately averages
fine state into coarse parents, combines both grids' maximum signal speeds
and diffusivities, and takes the minimum of their cut-cell hydro/diffusion
intervals, the user maximum and remaining time. The full interval is bounded
by the transport estimate, conservatively also bounding both half-steps.
Transactional rejection/halving remains available.

Fixed mode keeps V1 checkpoints. Automatic mode uses V2 magic with the same
record layout; each embedded patch config stores CFL controls. Restart rejects
mode or CFL mismatches. The retry limit remains a termination control rather
than numerical checkpoint context. Initial states are not changed by timestep
estimation. Empty fine patches contribute no bound.

These are pre-step estimates. They do not predict reaction-half-step changes,
control chemistry error, guarantee nonlinear coupled stability, or add
subcycling, multiple patches or MPI EB. The default remains fixed-cap behavior.

## Independently executed checks

GNU Fortran 13.3.0, Debug/Ninja, existing
`build/selected-tree-mpi-debug`; the exercised application is serial.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug \
  --target pelef_reactive_eb_two_level_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug \
  -R '^regression_reactive_eb_two_level_3d$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tools/check_reactive_eb_two_level_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final **1/1 passed**, zero failures/skips, 4.38 s (4.41 s total).
Provenance, Python syntax, project contract and whitespace checks passed.
Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

The checker retains existing application/restart/retry cases and adds:

- Coupled reaction/transport with automatic selection and a small user cap:
  final state and budget CSVs remain byte-identical to the fixed-cap reference.
- Hydro at final time 1e-4 s, user maximum 1e-3 s: 10 accepted steps, no
  rejected trials, normalized inventory residual `2.4368433926435087e-16`.
- First selected interval replayed directly in fixed mode: both CSVs match
  byte-for-byte. Restart after that step also exactly matches uninterrupted
  automatic execution. Changed mode or hydro CFL is rejected before output.
- Micrometre geometry with molecular transport: first interval
  `1.8758065462458827e-12` s. Halving only transport CFL halves the accepted
  interval, with no rejected trials in either run.

Logs: `build/eb-two-level-auto-dt-checks/{build,ctest,final-ctest}.log`.
The final log owns the counts; the earlier successful run preceded explicit
zero-rejection assertions for the diffusion cases. Artifacts:
`build/selected-tree-mpi-debug/tests/reactive_eb_two_level_3d/`, including
`auto.chk`, `auto_resume.nml.log`, `auto_wrong_cfl.nml.log`, and
`auto_diffusion_half.nml.log`.
No broader configuration matrix, installed audit, external solver parity or
physical validation was rerun.
