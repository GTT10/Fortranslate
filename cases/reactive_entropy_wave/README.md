# Reactive general-EOS entropy wave

This case transports a smooth density wave through a periodic domain at constant pressure, constant temperature, and pure N2 composition. Chemistry is disabled. It isolates the NASA7 conserved/primitive conversion, frozen-composition characteristic PLM, selectable Rusanov/PeleC-style fluxes, and the composition-dependent CFL path.

```bash
./build/pelef_reactive_1d \
  cases/reactive_entropy_wave/entropy_wave.nml
```

The automated regression runs both Rusanov and the qualified PeleC-style acoustic solver on 40, 80, and 160 cells and requires approximately second-order density convergence.
