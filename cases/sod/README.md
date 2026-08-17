# Sod shock tube

Canonical one-dimensional Sod Riemann problem:

| Region | Density | Velocity | Pressure |
|---|---:|---:|---:|
| Left, `x < 0.5` | 1.0 | 0.0 | 1.0 |
| Right, `x >= 0.5` | 0.125 | 0.0 | 0.1 |

The default run advances an ideal gas with `gamma = 1.4` to `t = 0.2` on 400 uniform cells. The solver writes `sod.csv` in the current working directory.

```bash
./build/release/pelef cases/sod/sod.nml
python3 tools/compare_sod.py --input sod.csv
```
