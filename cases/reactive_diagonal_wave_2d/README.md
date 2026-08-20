# Reactive diagonal waves in 2D

`diagonal_wave.nml` transports a periodic constant-pressure,
constant-composition entropy wave obliquely through the uniform Cartesian mesh.
The exact density is a translated sinusoid.

`diagonal_composition_ppm.nml` transports a constant-temperature,
constant-pressure H2/N2 composition wave with the characteristic-PPM normal
predictor and CTU correction. Density changes consistently with the local
mixture molecular weight. Together these cases verify density and species
transport independently of chemistry.
