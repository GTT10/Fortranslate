# Pinned Simple transport kernel — 2026-09-27

Implementation: `3d16d1b7`.

The existing analytical transport path uses Wilke viscosity and Mathur
conductivity mixing. Pinned PelePhysics `Source/Transport/Simple.H` instead
evaluates cubic log-temperature fits, a sixth-power viscosity mixture and a
fourth-root conductivity mixture. These are different algorithms; the
[evaporation discrepancy](evaporation-full-reference-20260927.md) cannot yet
be attributed quantitatively to this difference.

`pelephysics_simple_transport` now implements the ideal-gas Simple kernel in
Fortran, accepting molecular weights and upstream fit coefficients explicitly.
It retains the upstream mass-fraction trace regularization, triangular binary
fit ordering and species density-weighted diffusion formula. Outputs are SI
viscosity, conductivity, rho-D, and binary diffusion at the supplied pressure.
Invalid inputs return false and zero outputs. Existing transport selection and
all application call paths remain unchanged: fit import and application
connection are the next required step. No generated rates were edited.

## Completed checks

GNU Fortran 13.3.0 Debug rebuild of
`build/heptane-spray-mpi-debug/pelef_mpi_reactive_eb_tree_3d_selected_spray`
completed (log `build/simple-transport-build.log`). Configuration remains
MPI/spray ON, tests/SUNDIALS OFF, OpenMPI 4.1.6, selected heptane_3sp bundle.

The focused checker compiles and calls actual pinned C++ `SimpleTransport`
with heptane_3sp generated coefficients, exports those coefficients and calls
the Fortran routine. The reference is the existing GNU C++ 13.3.0 serial
Release spray build, PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8` and
PelePhysics `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-spray-mpi-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
python3 tools/check_simple_transport_reference.py --reference-build build/pelec-pinned-spray --candidate-build build/heptane-spray-mpi-debug --output build/simple-transport-reference-audit
```

Six states cover 300–1200 K, 101325 Pa/3 MPa, mixed gas and each pure species.
All 30 viscosity/conductivity/rho-D values pass the preset `2e-12` relative
gate; maximum relative error is `8.475e-16`. One negative-composition rejection
also passes. Binary diffusion positivity and symmetry pass; binary coefficient
magnitudes were not separately compared with a C++ API. Fortran checking uses
`-fcheck=all -ffpe-trap=invalid,zero,overflow`.

Artifacts: `build/simple-transport-reference-audit/` contains executable
harnesses, coefficients, states, outputs, exact commands and comparison JSON
with the source-export manifest and reference header hashes. The first
`build/simple-transport-reference/` run also passed; the second run adds
explicit source identity and finite-output checks to the checker.

Mechanism provenance, Python syntax, project contract and diff checks passed.
No CTest matrix or new full spray trajectory run was performed. This evidence
covers the transport kernel, not application integration or physical validity.
