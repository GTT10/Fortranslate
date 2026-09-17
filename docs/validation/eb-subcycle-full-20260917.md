# Coupled two-level EB application subcycling — 2026-09-17

Implementation: `4b9d36b1c8291f27aea685f635b344bfd882f542`, based on
`74c2c09` on local `main`.

The full and reaction/hydro library APIs accept optional `subcycle`. When true,
the composition uses ratio-based hydro and ratio-squared transport fine steps.
Reaction halves still solve composite active coarse/fine cells around the
complete `R/2 T/2 H T/2 R/2` sequence. Initial/accepted parent averaging, chemical
inventories, boundary/wall ledgers and limiter minima keep their existing
transactional contract. Disabled operators use the selected remaining operator
path. This symmetric arrangement still contains first-order Euler components.

The application exposes `subcycling` (default false). With automatic timestep
selection, a fine hydro interval permits at most `ratio` times that interval
on the coarse clock; a fine diffusion interval permits `ratio**2` times it.
Coarse bounds, user cap and final time still apply. Scaling uses capped
comparisons to avoid multiplying an unlimited interval. Both transport halves
remain conservatively bounded by the full-interval selection.

Checkpoint V1/V2 retain common-step fixed/automatic modes; V3/V4 identify the
corresponding subcycled modes with the same record layout. Restart requires
matching mode and existing physical/CFL context. This adds no multiple patches,
MPI EB, convergence proof or physical validation.

## Independently executed checks

GNU Fortran 13.3.0, Debug/Ninja, existing
`build/selected-tree-mpi-debug`; both tests are serial.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug \
  --target test_eb_two_level_full_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug \
  -R '^unit_eb_two_level_full_3d$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug \
  --target pelef_reactive_eb_two_level_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug \
  -R '^regression_reactive_eb_two_level_3d$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tools/check_reactive_eb_two_level_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Two focused invocations each passed **1/1**, zero failures/skips:
full composition 17.45 s (17.48 s total), application 6.33 s (6.35 s total).
Provenance, Python syntax, contract and whitespace checks passed.
Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

The full unit now runs common and subcycled modes with explicit/implicit
reaction, all-component boundary/wall/chemical balances, elemental conservation,
covered-parent masking, positive temperatures, disabled-operator exact parity,
and rollback including a hydro rejection after the first reaction and
transport halves.

The application checker adds a three-step coupled subcycled run and stopped
continuation. State and inventory CSVs are byte-identical after restart;
changing subcycling mode rejects before output. Normalized composite residual:
`4.1992326640465177e-16`, no rejected trials. Separate hydro/diffusion first-step
cases verify that fine-grid constraints permit larger coarse intervals bounded
by factors 2 and 4. The subcycled diffusion case accepts
`7.50322618498353073e-12` s, with reported minimum species limiter
`0.193995634123798416`. Fixed-cap V3 serialization is implemented but not
separately exercised by the new checkpoint regression; V4 continuation is.

Logs: `build/eb-subcycle-full-checks/{build,ctest,app-build,app-ctest}.log`.
Artifacts: `build/selected-tree-mpi-debug/tests/reactive_eb_two_level_3d/`,
including `subcycle.chk`, `subcycle_resume.nml.log` and
`subcycle_diffusion_stop.nml.log`.
No full matrix, installed audit, external solver parity, convergence study or
physical validation was rerun.
