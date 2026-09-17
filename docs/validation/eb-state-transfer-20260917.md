# Reactive EB coarse/fine state transfer — 2026-09-17

## Source and API

Implementation `82547ac13b2063031f337c522231e6232d47140d`, based on
`2d4daae` on local `main`.

`eb_reactive_transfer_3d_mod` provides full-domain
`restrict_reactive_eb_state_3d` and `prolong_reactive_eb_state_pcm_3d`.
Both take species, fine geometry and integer refinement ratio. The coarse
geometry is constructed using the conservative geometry restriction API, so
coarse fluid volume equals the sum of its fine fluid volumes.

Restriction averages conserved components using fine fluid-volume fractions.
Every active donor passes the existing conserved-to-primitive checks. The
volume-weighted recovered donor temperature supplies only an initial guess;
the destination temperature is recovered from its averaged conserved state.
PCM prolongation copies each coarse conserved state to its active fine children
and recovers temperature from the coarse state. It is first order in space and
preserves the same fluid-volume-weighted inventory.

Both operations preserve covered destination storage, including NaN payloads,
and ignore covered source storage. Candidate conserved fields and temperatures
are published together only after all donors and destinations succeed. Invalid
geometry, ratio, shape or active state leaves both destination arrays unchanged.
The geometry builder retains its existing allocation-exhaustion behavior.

This is a library transfer API. Patch subsets, limited higher-order interpolation,
coarse/fine ghost filling, EB reflux, hierarchy stepping and MPI integration
remain open. Callers must use the implied restricted coarse geometry; the API
does not accept an independently sampled coarse level set.

## Independently completed verification

GNU Fortran 13.3.0 Debug, Open MPI 4.1.6, Ninja; existing matching build
`build/selected-tree-mpi-debug`, MPI enabled, spray disabled. Core and selected
EB runtime compiled; execution used the fixed elementary mechanism. No MPI or
application matrix was run. No Python source changed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 2 --target test_eb_reactive_transfer_3d pelef_selected_reactive_eb_3d_runtime
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_eb_reactive_transfer_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 PASS, 0 failures, 0 skips**, test 0.06 s, CTest total 0.09 s.
The test uses a 6-cubed oblique plane and faceted sphere with spatially variable
density, pressure and velocity, each restricted/prolonged at ratios two and
three. It verifies conservation of every conserved component (including all
species), recovered temperature, coarse-to-fine-to-coarse consistency and bitwise
preservation of covered storage. Late invalid donors in each direction test
whole-operation rollback, together with invalid ratios and shape rejection.

Provenance and source contract passed. Mechanism SHA remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Logs: `build/eb-transfer-checks/{build,ctest}.log`.
These focused checks do not establish EB-AMR integration, convergence, external
PeleC parity or physical validation.
