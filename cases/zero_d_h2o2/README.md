# Elementary H2/O2 constant-volume reactor

This case uses seven species and four reversible elementary reactions selected from Cantera's `h2o2.yaml`. Reaction data are stored in `mechanisms/h2o2_elementary.json`, converted to SI units, and used to generate the committed Fortran mechanism module.

The initial mixture is stoichiometric H2/O2 diluted with N2 and contains a small radical seed. The case is an elementary-kinetics and Cantera-parity gate; it is not the complete Cantera H2/O2 mechanism because third-body and falloff reactions are intentionally excluded.

```bash
./build/pelef0d_h2o2 cases/zero_d_h2o2/reactor.nml
python3 tools/check_zero_d_h2o2.py --input zero_d_h2o2.csv
```

With Cantera 3.2 installed:

```bash
python3 tools/compare_h2o2_cantera.py \
  --input zero_d_h2o2.csv \
  --mechanism mechanisms/h2o2_elementary_cantera.yaml
```
