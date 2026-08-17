# Sod shock tube

Canonical one-dimensional Sod Riemann problem:

| Region | Density | Velocity | Pressure |
|---|---:|---:|---:|
| Left, `x < 0.5` | 1.0 | 0.0 | 1.0 |
| Right, `x >= 0.5` | 0.125 | 0.0 | 0.1 |

All supplied runs use `gamma = 1.4`, 400 uniform cells, CFL `0.45`, and final time `t = 0.2`.

## Piecewise-constant Rusanov baseline

```bash
./build/pelef cases/sod/sod.nml
python3 tools/compare_sod.py --input sod.csv
```

## PLM/MC with Rusanov

```bash
./build/pelef cases/sod/sod_plm.nml
python3 tools/compare_sod.py \
  --input sod_plm.csv \
  --density-l1-max 4e-3 \
  --pressure-l1-max 3e-3
```

## PLM/MC with the PeleC-style solver

```bash
./build/pelef cases/sod/sod_pelec.nml
python3 tools/compare_sod.py \
  --input sod_pelec.csv \
  --density-l1-max 2e-3 \
  --pressure-l1-max 1.5e-3
```

The verified GNU Fortran 14.2 Debug run produced density and pressure L1 errors of approximately `1.3678e-3` and `8.0552e-4`, respectively, while retaining positive states and roundoff-scale integral-balance errors.
