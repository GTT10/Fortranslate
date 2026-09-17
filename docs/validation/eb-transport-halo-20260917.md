# EB molecular-transport face donors — 2026-09-17

Implementation `ec5308392c7add67d8298f3ec99e9b2b98cfbf47`, based on
`48e8515` on local `main`.

## API and scope

`reactive_eb_transport_fluxes_rhs_3d` accepts optional conserved-state and
temperature face halos with spatial indices 0:n+1, plus a six-entry optional
`boundary_faces` mask in x-, x+, y-, y+, z-, z+ order. Both halo arrays must be
supplied together. With halos, all faces default enabled; without halos, all
default disabled. Enabling a face without halos is rejected.

On each enabled open boundary face, the donor is converted using the existing
EOS and passed to the existing directional transport kernel together with the
interior primitive state. Normal gradients use one ghost-center spacing;
tangential velocity gradients use the adjacent interior active-cell stencil.
Viscosity, conduction, species diffusion and barodiffusion flags are retained.
Boundary species flux and its enthalpy contribution enter the existing joint
species limiter before the aperture-weighted divergence is formed.

Disabled physical faces retain zero transport flux. Calls without the new
optional arguments retain the existing path. This is a boundary-flux interface,
not yet a synchronized hierarchy transport integrator. Its ghost-center/PCM
construction and one-sided tangential stencil do not establish higher-order
coarse/fine diffusion accuracy. Transport timestep selection, Euler/RK wrappers,
reflux of accepted transport flux and reacting hierarchy composition still need
to be connected. Existing failure-output semantics of the flux/RHS routine are
unchanged; callers must discard outputs when `ok` is false.

## Completed evidence

GNU Fortran 13.3.0 Debug, Open MPI 4.1.6, Ninja; existing matching
`build/selected-tree-mpi-debug`, MPI on, spray off. Fixed core and selected EB
runtime compiled. No MPI or application matrix was run. No Python source changed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 2 --target test_eb_transport_halo_3d test_reactive_eb_transport_3d pelef_selected_reactive_eb_3d_runtime
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^(unit_eb_transport_halo_3d|unit_reactive_eb_transport_3d)$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Final **2/2 PASS, 0 failures, 0 skips, total 1.17 s**:

- Existing coupled transport: 0.91 s.
- Transport halo: 0.22 s. Regular and cut-plane 4-cubed geometries exercise each
  enabled axis/side against independent directional face-kernel calls using
  different density, velocity, pressure, temperature and composition in the
  donor. Every conserved component satisfies the aperture-weighted boundary
  inventory identity. Large trial intervals activate the species limiter and
  leave nonnegative species inventories. Disabled faces retain the uniform
  zero-flux result; incomplete halo inputs are rejected.

The initial build failed because the new test used the size of an array in the
same ALLOCATE statement that created it. Splitting the allocation fixed it.
The premature initial CTest had one existing test pass and the new test Not Run;
it is not counted as successful verification. Final build and test logs are
`build/eb-transport-halo-checks/{final-build,final-ctest}.log`; initial failures
remain in `{build,ctest}.log`.

Provenance and source contract passed. Mechanism SHA remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
This focused result does not prove hierarchy transport integration, convergence,
external solver parity or physical validation.
