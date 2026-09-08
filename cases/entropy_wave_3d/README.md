# Three-dimensional entropy wave

This case advances a smooth density wave through a uniform periodic cube at
constant pressure and velocity. It is the public application case for the
first conservative three-dimensional PeleF update.

From a configured build directory, run:

```sh
./pelef3d ../cases/entropy_wave_3d/entropy_wave.nml
```

The application writes `entropy_wave_3d.csv` with x-fastest deterministic
cell ordering and reports the density L1 error and the change in every Euler
conserved integral. The `0.203.0` qualification uses the Rusanov flux, PCM
spatial states, SSPRK2 time integration, and periodic boundaries.
