# Zero-dimensional isomerization reactor

This is a deliberately small constant-volume reactor verification problem. Two equal-molecular-weight synthetic NASA7 species undergo the mass-conserving first-order conversion `A -> B`. Species B has a lower formation internal energy, so the adiabatic reactor heats while total specific internal energy remains fixed through the shared `e -> T` inversion.

```bash
./build/pelef0d cases/zero_d_isomerization/reactor.nml
python3 tools/check_zero_d_isomerization.py \
  --input zero_d_isomerization.csv
```

The default case starts at 700 K and approaches 1840 K as A is consumed. It is a thermodynamics and reactor-integration gate, not a detailed combustion mechanism or a claim of PeleC chemistry parity.
