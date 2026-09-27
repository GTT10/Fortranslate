# Fitted face transport coefficients — 2026-09-27

Implementation: `1ca9553c`.

Pinned PeleC `PeleCUtilities.H::pc_move_transcoefs_to_ec` defaults to a harmonic
mean of adjacent cell transport coefficients. The Fortran face helper instead
evaluated transport at the arithmetic mean temperature/composition/pressure.
For fitted Simple transport it now evaluates each cell independently and
harmonically averages viscosity, conductivity and density-weighted diffusion.
The returned diffusion coefficient is divided by face density to preserve the
existing flux API. The analytical backend retains its previous behavior.

This shared face helper serves 2D and 3D fluxes; bulk viscosity remains connected
only in the previously qualified 3D path. Fitted checkpoint identity advances
to `PELEF_SIMPLE_TRANSPORT_V3`, distinguishing the changed face model. No
coefficient values or comparison tolerances were altered.

## Completed checks

Existing GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6 configuration,
MPI/spray ON, tests/SUNDIALS OFF:

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
python3 tools/check_simple_transport_reference.py --reference-build build/pelec-pinned-spray --candidate-build build/heptane-simple-mpi-debug --bundle build/heptane-3sp-simple-reference/heptane_3sp.json --bulk --output build/simple-face-record-reference
```

The final build passes (`build/simple-face-final-build.log`). The initial build
failed at linking because the new calls used a helper name from the 3D module;
they were corrected to the existing 2D `multiply_finite` helper. That failure
is retained in `build/simple-face-build.log` and is not counted as a pass.

The compiled reference check covers six states / 36 cell coefficients and three
nonuniform faces / 15 face coefficients. Reference face expectations use the
pinned harmonic formula applied to independently compiled C++ cell coefficients.
Maximum face relative error: `6.159e-16`; the unchanged gate is `2e-12`.
Cell coefficients remain within `1.103e-15`. The compression-stress, viscous-work
and timestep checks continue to pass. Fitted checkpoint write/read round-trip
and rejection of the preceding V2 marker pass. This is one focused compiled
harness, not a CTest matrix. Exact compiler commands and hashes are in
`build/simple-face-record-reference/`; the preceding face-only extension check
is `build/simple-face-reference/`.

One two-step evaporation application run also completes:

```sh
python3 tools/prepare_evaporation_reference.py --bundle build/heptane-3sp-simple-reference/heptane_3sp.json --output build/spray-face-transport
# From the prepared directory:
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none /usr/bin/mpiexec -n 1 /home/gtt13/projects/pelec/build/heptane-simple-mpi-debug/pelef_mpi_reactive_eb_tree_3d_selected_spray tree.nml spray.nml
```

Against the unchanged full C++ evaporation reference, axial momentum still
differs by `2.654e-6`, total energy by `2.360e-10`, and stored temperature by
`2.218e-9`, normalized by each reference field maximum. Fuel-density difference
is `1.975e-12`. The face change therefore does **not** explain this case's
remaining momentum/temperature residual. Artifact:
`build/spray-face-transport/full-comparison.json`, including input hashes.

Provenance, Python syntax, project contract and diff checks pass. No new
rank-changing application restart, installed-binary audit, full test matrix or
physical validation was run. Full spray parity remains open.
