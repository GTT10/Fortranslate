# Pinned linear momentum forcing — 2026-09-27

Implementation: `3206d77ba1a2ada5106732e836f148336803993a`.
Regular fixed/selected SDC applications now accept
`reactive_linear_forcing` with `coefficient` and three `momentum_offset`
values. Defaults are zero. Nonzero coefficient requires SDC.

The source follows pinned PeleC `Source/Forcing.cpp` at
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8` literally:

`momentum_rate[d] = coefficient * rho * (momentum[d] - momentum_offset[d])`.

Density/species/total-energy rates are zero. In particular, the subtracted
quantity is momentum density, despite upstream's `forcing_u0/v0/w0` names;
this is not a velocity relaxation law. No mechanical-work source is added.
The ordinary external-force/work model remains a distinct source and can
coexist. SDC evaluates the state-dependent forcing at each source stage and
includes it in the source inventory. Transport and LES can participate too.

PeleF uses SI: coefficient m³/(kg s), offsets kg/(m² s). To reproduce upstream
CGS numbers, coefficient_SI = 0.001 coefficient_CGS and offset_SI = 10 offset_CGS.
Upstream HIT also uses its named offsets as initial mean velocities; that
separate initialization must be configured explicitly in PeleF for nonzero
offset comparisons. The checks here use zero offsets.

Parameters must be finite with magnitude at most `1e30`; source evaluation
checks density/finite bounded momentum before multiplication. Failed candidates
retain the existing SDC transaction. Other application families do not yet
expose this input.

## Completed checks

Unchanged GNU Fortran 13.3.0 Debug MPI/spray/SUNDIALS configuration; rebuilt
`pelef_reactive_3d` successfully (`build/hit-linear-forcing-build.log`).
Pinned full reference/exporter identities remain those in
[the reference record](hit-full-reference-20260927.md).

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/amr-cvode-app-debug --target pelef_reactive_3d --parallel 2
python3 tools/check_hit_evolution_reference.py --pelec build/pelec-pinned-hit/Exec/RegTests/HIT/PeleC-HIT --fortran build/amr-cvode-app-debug/pelef_reactive_3d --exporter build/pelec-plotfile-export/pelef_plotfile_csv --velocity build/pelec-pinned-hit-source/Exec/RegTests/HIT/hit_ic_4_32.dat --work build/hit-full-reference --steps 3 --linear-forcing 1000
```

Repeated with `--transport --les wale`. Two successful comparisons, 4³ cells,
three `1e-8` s steps each; zero failures/skips. Both use coefficient 1000 SI
(1e6 CGS). Maximum normalized field difference is `3.916735033144481e-15`,
below the unchanged `1e-10` gate. Source-accounted Euler residuals are
`1.1288783166368997e-16` and `2.2577566391289868e-16`, respectively.
Momentum source inventories are nonzero (approximately `[-1.61,-2.74,1.68]e-8`);
the hydro-only energy source inventory is exactly zero.

Artifacts under `build/hit-full-reference/`:

- Hydro + forcing: `evolution_54cbie_s/`.
- Transport + WALE + forcing: `evolution_uswmrtqx/`.

Each retains logs, inputs, plotfiles, CSVs, field errors and binary hashes.
A separate `coefficient=1e31` input was rejected before CSV publication
(`build/hit-linear-forcing-reject/`). Mechanism provenance, Python syntax,
project contract and diff checks passed. No CTest matrix, selected executable,
nonzero-offset case or long-time forced-turbulence statistics were rerun.
