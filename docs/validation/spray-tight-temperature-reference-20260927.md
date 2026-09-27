# Full evaporation reference with tighter temperature inversion — 2026-09-27

This follows the [compiled inversion audit](spray-temperature-inversion-20260927.md).
It is a diagnostic variant of the pinned reference, not a replacement baseline.
No Fortran numerical implementation or comparison tolerance changes.

Pinned PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`, PelePhysics
`6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`, AMReX
`1546c1ccbf41c704a33790f4bf445a2fcee1372e`, SUNDIALS
`0eff39663606f2ff280c4059a947ed62ae38180a`. The Fortran output remains the
`595bc61e` harmonic-face run in `build/spray-face-transport/`; the subsequent
default-control check established byte-identical output at `46db17e1`.

## Completed intervention

GNU C++ 13.3.0 Release, serial CPU, particles ON, tests OFF, fresh directory:

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S build/pelec-pinned-hit-source -B build/pelec-pinned-spray-tight -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DCMAKE_CXX_FLAGS=-DCONVERGENCE -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DPELE_ENABLE_MPI=OFF -DPELE_ENABLE_PARTICLES=ON -DBUILD_TESTING=OFF
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/pelec-pinned-spray-tight --target PeleC-Spray-Evaporation --parallel 2
# From build/spray-evaporation-tight-reference, with byte-identical reference.inp and particles.dat:
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ../pelec-pinned-spray-tight/Exec/RegTests/Spray-Evaporation/PeleC-Spray-Evaporation reference.inp
```

The source's existing `CONVERGENCE` branch changes temperature inversion from
200 iterations / `1e-6 K` stopping correction to 5000 / `1e-12 K`.
The full build and the same 4³ two-step run complete. Logs:
`build/pelec-pinned-spray-tight-{configure,build}.log` and
`build/spray-evaporation-tight-reference/reference.log`.
The run directory records input hashes, binary SHA-256, exported CSVs and
`comparison.json`. Initial gas CSVs are byte-identical between references.

| Difference from unchanged Fortran output | Default reference | Tight reference |
| --- | ---: | ---: |
| Maximum temperature, K | `8.87208e-7` | `3.75127e-9` |
| Axial momentum / reference field maximum | `2.65361e-6` | `2.65863e-6` |
| Total energy / reference field maximum | `2.35965e-10` | `2.34892e-10` |
| Fuel density / reference field maximum | `1.97462e-12` | `1.97462e-12` |

The stopping rule explains most of the stored-temperature residual but does
not explain the remaining momentum difference. This is an intervention result,
not a full parity acceptance test or mesh/time convergence study.

## Energy accounting

Using exported initial/final gas energy and particle diameter, temperature,
velocity and multiplicity, the combined reference energy changes by
`-2.50951e-15 J` (default) and `-2.50946e-15 J` (tight). The corresponding
Fortran application budget residual is `1.89735e-19 J`.

The reconstruction uses cell volume `1e-9 m3`, liquid density `681.4 kg/m3`,
and liquid specific energy
`h_NC7H16(300 K) - 363000 + 2248.3*(T-300) J/kg`, plus kinetic energy.
Cantera 3.2.0 supplies the reference enthalpy from the pinned YAML. The saved
script, input hashes and results are in
`build/spray-temperature-kernel-audit/{check_energy.py,reference-energy-audit.json}`:

```sh
build/cantera-venv/bin/python build/spray-temperature-kernel-audit/check_energy.py
```

Pinned `Drag.H` accumulates explicit drag/convection sources and final-state
vapor enthalpy. Fortran's mixed-particle routine instead returns the negative
of the actual particle energy change. The measured accounting difference is
consistent with that discretization distinction; it does not by itself prove
the spatial momentum residual's cause. Retain conservative Fortran feedback.

Provenance, Python syntax, project contract and diff checks pass. No broad
matrix, new candidate simulation, physical validation or acceptance-gate
change is included.
