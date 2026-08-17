# Shu-Osher shock-density-wave interaction

This case uses the canonical one-dimensional Shu-Osher problem on
`[-5, 5]`. A Mach-like shock state begins to the left of `x = -4`,
while the right state contains the sinusoidal density field

\[
\rho(x)=1+0.2\sin(5x), \qquad u=0, \qquad p=1.
\]

The supplied configuration uses 800 cells, PLM with the MC limiter, the
single-species PeleC-style approximate Riemann solver, CFL `0.25`, and
final time `t = 1.8`.

```bash
./build/pelef cases/shu_osher/shu_osher.nml
python3 tools/check_shu_osher.py --input shu_osher.csv
```

The regression checks positivity, integral balances including boundary
fluxes, retained oscillatory structure, and two pinned density-field
signatures. It is a deterministic regression, not an exact-solution
comparison.
