# EB MOL pressure/flux diagnostic output — 2026-10-23

## Scope

This is diagnostic infrastructure for the outstanding transverse shear maxima,
not a numerical correction or a C2 qualification result. The next step is to
capture these outputs at identified recursive stages and compare the x/cut-plane
and x/z coarse/fine intersections before and after tree face synchronization.
No new shear campaign has been launched with this interface yet.

`mol_riemann_flux_3d` has optional `reconstructed_pressure(2)` output containing
the actual left/right pressures consumed by the solver, after characteristic
reconstruction and positivity limiting. Failed calls preserve the output.

`compute_mol_eb_hydro_3d` has optional `eb_mol_diagnostics_3d` output:

- `cell_pressure` is the EOS primitive pressure on the complete input box.
- `pressure` holds the two reconstructed pressures at each Cartesian face.
- `raw_flux` holds the flux before EB centroid interpolation and limiting.
- `centroid_flux` holds the flux used in the local divergence, after limiting.
- `rate` and `wall_rate` hold the local spatial rate and wall contribution.

Raw and centroid flux arrays include an extra final internal-energy component
when the auxiliary energy equation is requested; conserved components are 1:nv.
All fluxes in the diagnostic are instantaneous, unlike the ordinary
`integrated_flux` return. The face-shaped container is reused only for layout.
Raw faces and reconstructed pressures include the transverse one-cell margin:
face index zero corresponds to input index 2 in its normal direction; other
indices map to input index +2. Centroid faces and rates use the physical
box at input cells 4:n-3. With six input halos this includes the spatial-rate
halo, not only the valid patch. Pressure at closed faces is zero and must be
masked by geometry; it is not a physical zero-pressure state.

The optional snapshot is published only after successful completion. No global
trace state or file I/O is introduced. Stage/node/time identity and output
selection belong to callers. Tree-level boundary limiting, synchronization,
reflux and accepted time increments are not represented by this local snapshot.

## Verification status

Base: dirty `main`, HEAD `0a0e00ebe9138145fbb0e0e745665d0d768df009`;
HEAD tree `f67d918cf8d0ee3a3c20d330c369a844cb591860`. Pre-edit copies of the
four changed files are in `build/eb-stage-flux-20261023/`; existing user changes
were preserved. This work does not identify the dirty source by HEAD alone.

Build uses the existing GNU Fortran 13.3.0 Debug/OpenMPI configuration at
`build/eb-eos-application-20261021/isolated-openmpi/registration`.
Targets: `test_eb_mol_hydro_3d test_mol_riemann_3d`, parallelism 4.
Build log: `build/eb-stage-flux-20261023/build.log`.

Tests extend the existing fixtures to check actual reconstructed pressure,
unchanged operator outputs, instantaneous flux units, directional momentum
closure including cut-cell wall contributions, halo indexing and rollback.
Completed focused CTest: **2/2 passed**, zero skips, 3.14 seconds:
`ctest --test-dir build/eb-eos-application-20261021/isolated-openmpi/registration
-R '^unit_(eb_mol_hydro_3d|mol_riemann_3d)$' --output-on-failure`.
JUnit and final log are `tests.xml` and `tests.log` in the artifact directory.
The initial run failed the newly added halo test because its comparison included
the auxiliary energy component against an nv-component conserved array. The
comparison now explicitly selects 1:nv; the failed log is retained as
`tests-first-failed.log`. No numerical tolerance was relaxed.
Mechanism provenance, project contract, Python syntax and diff whitespace
checks passed after the final code edits. No application/MPI run or fixed-PeleC
comparison was performed in this increment. Stage capture remains unfinished.
