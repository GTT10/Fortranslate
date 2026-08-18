# Passive two-species Sod problem

This case uses the existing constant-`gamma` ideal-gas hydrodynamics and adds two
conserved species densities. Species 1 labels the initial left state and species
2 labels the initial right state. The milestone verifies conservative species
transport and mass-fraction closure; it does not yet model species-dependent
thermodynamics.

```bash
./build/pelef_ms cases/multispec_sod/multispec_sod.nml
python3 tools/check_multispec_sod.py --input multispec_sod.csv
```
