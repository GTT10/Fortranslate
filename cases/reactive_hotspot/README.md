# Reactive H2/O2 hotspot

A Gaussian temperature hotspot is initialized at uniform pressure in a
periodic stoichiometric H2/O2/N2 mixture with a small radical seed. The case
exercises NASA7 state recovery, frozen-composition characteristic PLM,
selectable general-EOS Riemann fluxes, and Strang-split cell chemistry.

The robustness baseline is:

```bash
./build/pelef_reactive_1d cases/reactive_hotspot/hotspot.nml
```

The qualified PeleC-style acoustic-flux regression is:

```bash
./build/pelef_reactive_1d cases/reactive_hotspot/hotspot_pelec.nml
```
