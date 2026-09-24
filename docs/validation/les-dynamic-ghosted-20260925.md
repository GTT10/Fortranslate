# Dynamic LES supplied-halo patch evaluation — 2026-09-25

Implementation `8cc44d26e0368b83b361bbabb88fa96d82caa771`, based on `d056d2b`.

The new `compute_dynamic_les_rhs_ghosted_3d` API accepts caller-owned donor halos
and returns an interior RHS and all three face-flux arrays transactionally.
Its dependency radius is test-filter radius plus five: derivative radius one,
coefficient smoothing radius three and adjacent face cell radius one. It rejects
insufficient halos, inconsistent dimensions and invalid states. The implementation
reuses the existing padded periodic operator; its seam is outside the interior
dependency stencil. All padded/filter states must be thermodynamically admissible.

This establishes patch-local input/output, not AMR Dynamic LES completion.
AMR donor construction, face synchronization/reflux, timestep selection and
coupled/application dispatch remain. The implementation also overcomputes padded
cells and does not yet expose an interior-only diffusion bound.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug, existing fixed tests-enabled and selected H2O2
tests-disabled configurations. Fixed `test_spray_core` and selected
`pelef_selected_spray_runtime` built; selected execution was not run.

One direct `test_spray_core les_dynamic_ghosted` invocation passed with exit 0.
The 4-cubed nonuniform periodic reference is embedded in an 18-cubed workspace
with seven donor layers. The outer extent is deliberately not divisible by the
period. Interior RHS and all faces match the reference exactly. Insufficient
halo and invalid interior gas state preserve all output arrays. Existing dynamic
flux conservation, translation, face-divergence and rollback checks also execute.
No tolerances changed.

Provenance, Python syntax, project contract and whitespace passed. The previous
MPI launch remains unresolved: CTest 2018710 / mpiexec 2018713 were still present
at this turn's final observation. No additional MPI launch or CTest run was
started. This result is a direct serial test, not MPI/AMR qualification.
No external solver parity, physical validation or broad matrix was run.

## Commands and artifacts

Logs: `build/les-dynamic-ghosted-checks/` (`provenance.log`, `build.log`,
`build-selected.log`, `test.log`, `contract.log`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_core --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH build/spray-amr-reactive-debug/tests/test_spray_core les_dynamic_ghosted
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
