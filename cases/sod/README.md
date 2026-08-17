# Sod shock tube

Canonical one-dimensional Sod Riemann problem:

| Region | Density | Velocity | Pressure |
|---|---:|---:|---:|
| Left, `x < 0.5` | 1.0 | 0.0 | 1.0 |
| Right, `x >= 0.5` | 0.125 | 0.0 | 0.1 |

Both supplied runs use `gamma = 1.4`, 400 uniform cells, CFL `0.45`, and final time `t = 0.2`.

## Piecewise-constant baseline

```bash
./build/pelef cases/sod/sod.nml
python3 tools/compare_sod.py --input sod.csv
```

This preserves the original first-order Rusanov regression baseline.

## Piecewise-linear reconstruction

```bash
./build/pelef cases/sod/sod_plm.nml
python3 tools/compare_sod.py \
  --input sod_plm.csv \
  --density-l1-max 4e-3 \
  --pressure-l1-max 3e-3
```

The PLM case reconstructs primitive variables with the monotonized-central limiter. In the verified GNU Fortran 14.2 Debug run, the density and pressure L1 errors were approximately `1.891e-3` and `1.198e-3`, respectively, while maintaining positive density and pressure and conservative integral balances.
