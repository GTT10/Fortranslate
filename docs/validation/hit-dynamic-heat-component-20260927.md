# Dynamic HIT heat-component difference — 2026-09-27

Implementation: `df012a3cd9400b6deba53f863556dedac07475e5`.
This resolves the cause of the bounded residual recorded in
[the preceding comparison](hit-dynamic-reference-20260927.md).

## Source finding and intervention

In pinned PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`,
`Source/LES.cpp` lines 579–583 pass component **0** of the cell-centered
three-component `flux_T` to every directional face interpolation.
`Source/LES.H` then uses this scalar in the face energy flux. PeleF's default
dynamic model uses the face-normal component instead.

The new optional `pinned_heat_flux=.true.` in `reactive_les` reproduces the
fixed reference's x-component selection for all directions. It requires
`model='dynamic'`. Default is false: the face-normal formulation remains the
normal model behavior. This option is for explicit reference compatibility,
not a claim that selecting x for every face is physically preferable.
Periodic, boundary and ghosted dynamic RHS entries forward the optional flag;
wall heat closures remain applied afterward. Only periodic regular application
behavior was exercised here. Other application wrappers do not expose it.

With the same initial field, transport, filters, timestep and binary on both
PeleF runs, toggling only this flag reduces the maximum normalized difference
from `7.67090568842832e-11` to `3.916735033144481e-15`. Temperature difference
falls to `9.473903143468002e-16`. This intervention supports the component
selection as the cause of the earlier difference for this case.

## Completed checks

Unchanged GNU Fortran 13.3.0 Debug MPI/spray/SUNDIALS configuration; rebuilt
`pelef_reactive_3d` successfully. Build log: `build/hit-pinned-heat-build.log`.
Reference/exporter binaries and pinned source identities remain those in
[the full-reference record](hit-full-reference-20260927.md).

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/amr-cvode-app-debug --target pelef_reactive_3d --parallel 2
python3 tools/check_hit_evolution_reference.py --pelec build/pelec-pinned-hit/Exec/RegTests/HIT/PeleC-HIT --fortran build/amr-cvode-app-debug/pelef_reactive_3d --exporter build/pelec-plotfile-export/pelef_plotfile_csv --velocity build/pelec-pinned-hit-source/Exec/RegTests/HIT/hit_ic_4_32.dat --work build/hit-full-reference --transport --les dynamic --pinned-heat-flux
```

Repeated without `--pinned-heat-flux` to verify unchanged default behavior.
Both comparisons passed the unchanged `1e-10` gate: two successful runs,
zero failures/skips. The default Fortran CSV was byte-identical to the
preceding `evolution_9m4m8ugc` result. Preflight, Python syntax, project contract
and `git diff --check` passed; no CTest matrix was rerun.

Artifacts under `build/hit-full-reference/`:

- Pinned selection: `evolution_myw6bzt6/`.
- Default normal selection: `evolution_tanlord3/`.

Each retains inputs, logs, plotfiles, CSVs and summaries with executable hashes.
Scope is still 4³ cells, one `1e-8` s step; no multiple-step or physical
turbulence validation is claimed.
