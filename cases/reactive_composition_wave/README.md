# Reactive general-EOS composition wave

This non-reacting periodic case advects a smooth H2/N2 composition wave at
constant pressure and temperature. The changing mixture molecular weight
produces a coupled density wave. It exercises the frozen-composition
characteristic PLM path and the selectable general-EOS HLLC flux.

```bash
./build/pelef_reactive_1d \
  cases/reactive_composition_wave/composition_wave.nml
python3 tools/check_reactive_composition_wave.py \
  --input reactive_composition_wave.csv
```

The regression also transports a discontinuous material interface and requires
HLLC to resolve it more sharply than the Rusanov baseline.
