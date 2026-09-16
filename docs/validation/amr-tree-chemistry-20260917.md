# Regular 3D tree chemistry and reaction/hydro split — 2026-09-17

Source parent: `57ca6b9d81ab798838d62ff6c344a0cd7e547ba3` on `main`, plus this change.
PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
New chemistry module SHA256:
`77f931dcc30430b49e89779456fb75581d9637b4b05dc32e5ad01110823e0166`.
New test SHA256:
`4a8312d52cdb7a0f0c7352be9c465bc96da7b0a964f97daa5b16d6fbede71546`.

## Implemented scope

The source API advances exposed cells with the existing chemistry solver,
restricts children bottom-up and recovers synchronized temperatures. Covered
parent cells do not repeat the chemistry solve. The coupled API composes
R/2-H-R/2 with synchronous tree hydro. Whole-hierarchy candidate copies preserve
state, temperature and flow diagnostics if any source/hydro/recovery phase fails.
Species/reactions and optional integrator are supplied by the caller; no fixed
database imports were added to the implementation.

This is a library path in `pelef_core`, exercised with elementary H2/O2.
It does not establish arbitrary-mechanism qualification, a selected-runtime
frontend, tree MPI/application/restart, molecular transport, subcycled time
integration, second-order coarse/fine spatial interpolation or physical validation.

## Focused verification

Reused `build/spray-mpi-migration-gnu-debug`: GNU Fortran 13.3 Debug with bounds
checks and invalid/zero/overflow traps, OpenMPI 4.1.6; MPI/spray enabled.
No toolchain change.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug \
  --target test_amr_patch_tree_chemistry_3d test_amr_patch_tree_hydro_3d -j 4
ctest --test-dir build/spray-mpi-migration-gnu-debug \
  -R '^unit_amr_patch_tree_(chemistry|hydro)_3d$' --output-on-failure
# After correcting the new test fixture, rebuild/rerun only that target:
cmake --build build/spray-mpi-migration-gnu-debug --target test_amr_patch_tree_chemistry_3d -j 4
ctest --test-dir build/spray-mpi-migration-gnu-debug \
  -R '^unit_amr_patch_tree_chemistry_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Hydro passed (8.01 s). The initial new chemistry test failed because its
radical-free elementary mixture had no measurable source activity at the chosen
interval. Adding a trace H fraction made the source active. A subsequent test
failure exposed a comparison against unaveraged parent initialization; the
oracle and initial hierarchy now follow the same existing restriction rule.
No production-code change or tolerance relaxation was needed for these failures.
The corrected chemistry test passed (33.24 s). Final evidence therefore covers
**two distinct passing tests**, not an uninterrupted 2/2 CTest invocation.

The new test uses three levels/four nodes, adjacent root children and a deeper
leaf. Both explicit and implicit chemistry are compared with the existing
uniform single-grid source result. A perturbed fine cell exercises nonuniform
periodic/outflow reaction-hydro evolution. Checks cover Euler invariants,
composite H/O/N balance including boundary flux, species closure, exact final
restriction consistency, invalid-reaction rollback and hydro failure after a
successful first half reaction step.

Provenance and project contract passed; both checker scripts passed Python AST
parsing. Whitespace checks passed. Logs are in `build/amr-tree-chemistry-checks/`:
`build.log`, `ctest.log` (initial invocation) and `chemistry-final.log`.
