# Full H2/O2 constant-volume reactor

This case exercises the ten-species, 29-reaction H2/O2 mechanism normalized from Cantera `h2o2.yaml`. It includes third-body reactions, duplicate reactions, and the Troe falloff reaction `2 OH (+M) <=> H2O2 (+M)`.

```bash
./build/pelef0d_h2o2_full cases/zero_d_h2o2_full/reactor.nml
python3 tools/check_zero_d_h2o2_full.py --input zero_d_h2o2_full.csv
```
