# Full pinned HIT executable comparison — 2026-09-27

The full pinned PeleC HIT executable and PeleF agree for one `1e-8` s step on
a periodic 4³ grid using the pinned 32³ velocity input. Three cases passed:
hydro only, hydro plus constant transport, and hydro/transport/WALE. This is
bounded evolved-field parity, not a long-time turbulence or physical validation.

## Source and build identity

- Fortran numerical implementation: `90085c9b2908540b24a4ca25e7981b9fe6e2e2ea`.
- Checker/exporter: `b52163c7`, extended by
  `efdc35c2f4f2f464dfd705d1b112a7d544921ed3` for transport/LES.
- PeleC: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
- PelePhysics: `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`.
- AMReX: `1546c1ccbf41c704a33790f4bf445a2fcee1372e`.
- SUNDIALS: `0eff39663606f2ff280c4059a947ed62ae38180a`.

Reference sources were exported using local `git archive` at each pinned commit,
without modifying checkouts or numerical source. The fresh build used GNU C/C++
13.3.0 Release, double precision, CPU serial, particles off, tests off. The
`PeleC-HIT` target completed successfully. Embedded Git labels are `unknown`;
the exported provenance manifest owns identity. Fortran reused the unchanged
GNU Fortran 13.3.0 Debug MPI/spray/SUNDIALS configuration. The serial CSV exporter
was built against this exact reference AMReX library.

Build commands are in the [case guide](../../cases/reactive_hit_3d/README.md).
Logs: `build/pelec-pinned-hit-configure.log`,
`build/pelec-pinned-hit-build.log`, `build/pelec-plotfile-export-configure.log`,
`build/pelec-plotfile-export-build.log`. No build failure occurred.

## Completed comparison

```sh
python3 tools/check_hit_evolution_reference.py --pelec build/pelec-pinned-hit/Exec/RegTests/HIT/PeleC-HIT --fortran build/amr-cvode-app-debug/pelef_reactive_3d --exporter build/pelec-plotfile-export/pelef_plotfile_csv --velocity build/pelec-pinned-hit-source/Exec/RegTests/HIT/hit_ic_4_32.dat --work build/hit-full-reference
```

Repeated once with `--transport` and once with `--transport --les wale`.
Three comparisons passed, zero failures/skips. Each compares 64 cells in
density, three momenta, total energy density, temperature and pressure, with
matching coordinates/time. PPM, artificial viscosity 0.1, two SDC iterations,
gamma 1.4, 300 K, 101300 Pa, Mach 0.1, 2π cm domain. Constant transport derives
from Taylor length 0.5 cm, Re 100, Pr 0.71. WALE uses Cw 0.5, CI 0.09, PrT 0.7.
Chemistry/species diffusion are disabled. CGS output is converted to SI.

| Case | Maximum normalized field difference | Maximum reference field update | Artifact directory under `build/hit-full-reference/` |
| --- | ---: | ---: | --- |
| Hydro | 3.916735033144481e-15 | 7.645655122947082e-4 | `evolution_cnbj_6l_` |
| Constant transport | 3.766091378023539e-15 | 7.653281857626366e-4 | `evolution_tmw2p2ak` |
| Transport + WALE | 3.916735033144481e-15 | 7.859034312064012e-4 | `evolution_e1uv8eoq` |

The predefined field gate is `1e-10`. Normalizers are initial density,
rho times target RMS velocity, initial internal-energy density, temperature,
and pressure. The reference must materially evolve (`max update > 1e-7`).
Artifacts retain both executables' logs, inputs, reference plotfiles, exported
CSVs, and summaries with binary and velocity-file SHA-256 hashes. No tolerance
or solver changes were needed to pass.

Mechanism provenance, Python syntax, project contract and `git diff --check`
also passed. No full CTest/configuration matrix was rerun. Smagorinsky/Vreman
options in the checker were not run here; dynamic LES, multiple timesteps,
larger grids, forcing, MPI/AMR and long-time statistics remain outside this
comparison.
