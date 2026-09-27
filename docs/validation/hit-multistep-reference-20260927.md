# Bounded timesteps and multistep HIT reference — 2026-09-27

Implementation: `5354d5377d2220c98c19542cc4feb5eb6e56d7e1`.
The regular fixed/selected application accepts optional
`&reactive_timestep maximum_dt=... /`. The finite positive limit caps the
already computed hydro, inlet, transport, wall and LES timestep. It cannot
override any stability bound. Omission adds no limit. A step that cannot
advance floating-point time is rejected. This input is specific to regular
applications; other front ends do not expose it.

The reference checker now accepts `--steps N` (1–10000), uses `1e-8` s as the
upper bound, and requires the requested completed step count and matching
final time. If stability requires smaller steps, the run fails rather than
silently comparing different trajectories. The final reference plotfile is
selected by step number. Existing field gate `1e-10` is unchanged.

## Completed checks

GNU Fortran 13.3.0 Debug, unchanged MPI/spray/SUNDIALS configuration. Rebuilt
`pelef_reactive_3d` successfully (`build/hit-multistep-build.log`). Pinned
PeleC/exporter identities and setup remain those in
[the full-reference record](hit-full-reference-20260927.md).

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/amr-cvode-app-debug --target pelef_reactive_3d --parallel 2
python3 tools/check_hit_evolution_reference.py --pelec build/pelec-pinned-hit/Exec/RegTests/HIT/PeleC-HIT --fortran build/amr-cvode-app-debug/pelef_reactive_3d --exporter build/pelec-plotfile-export/pelef_plotfile_csv --velocity build/pelec-pinned-hit-source/Exec/RegTests/HIT/hit_ic_4_32.dat --work build/hit-full-reference --transport --les dynamic --pinned-heat-flux --steps 3
```

Repeated for `--les smagorinsky` and `--les vreman`, omitting
`--pinned-heat-flux`. Three successful comparisons, each 4³ cells and three
steps, final time approximately `3e-8` s. Zero comparison failures/skips.

| Model | Maximum normalized field difference | Artifact directory under `build/hit-full-reference/` |
| --- | ---: | --- |
| Dynamic, pinned heat component | 3.916735033144481e-15 | `evolution_jrdajmxg` |
| Smagorinsky | 7.711757158782954e-14 | `evolution_c3mdlvj3` |
| Vreman | 3.012701199622825e-14 | `evolution_7u4f9s0d` |

Maximum reference field updates range from `2.34e-3` to `3.01e-3` in the
same normalization, demonstrating active evolution. Summaries record each
binary hash, field errors and updates; inputs/logs/plotfiles/CSVs are retained.
Smagorinsky/Vreman maxima are in temperature; momentum and total-energy
differences are below `4.4e-15` and `1.2e-15`, respectively.

A separate zero-limit input was rejected with `Invalid maximum timestep`
before CSV output (`build/hit-timestep-rejection/`). Provenance preflight,
Python syntax, project contract and diff checks passed. No CTest matrix or
selected executable was rerun. WALE's earlier one-step comparison remains
the evidence for that model. Longer-time statistics and physical turbulence
validation are outside these checks.
