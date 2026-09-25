# Single-level EB inlet compositions — 2026-09-26

Implementation 65a0827, parent d0d2394, local main, version 0.246.0.
Fixed and selected single-level gas applications now resolve named per-face
mole fractions through the active mechanism, convert to mass fractions and
construct inlet conserved states from density, pressure and velocity. Each
unspecified face independently retains initial gas composition.

Configuration validation rejects invalid fractions, duplicate/blank names,
unused data and overrides on non-inflow faces. The application additionally
rejects names absent from its mechanism. The optional INLET_COMPOSITION_V1
checkpoint section binds exact counts, species order and input fractions;
legacy inputs retain their prior checkpoint representation. Read failure uses
the existing candidate-state transaction. Single-level spray input rejects
these overrides until its state construction is connected.

## Completed checks

GNU Fortran 13.3.0, Python 3.12.3, Ninja Debug and existing OpenMPI wrappers.
Fixed tests-enabled and selected tests-disabled builds remained separate;
selected bundle is H2O2 full. No compiler/backend configuration changes.

- Existing configuration/checkpoint tests: **2/2 passed**, 0.71 s reported total.
- Fixed application: **6 successful executions, 10 expected rejections**.
- Selected application: **6 successful executions, 10 expected rejections**.
- Inert, reacting, coupled, checkpoint-stop, resumed and altered-mixture runs
  completed. Continuous/resumed CSVs are byte-identical; altered-mixture CSVs
  differ. Existing mass-fraction and component inventory checks pass.
- Rejections cover method/solver configuration, invalid density, changed inlet
  velocity and labels, changed/removed restart composition, unknown species and
  negative fractions. Rejected processes did not produce their requested CSVs.

First checker attempts on both mechanisms stopped at the removed-composition
assertion: the fixture reset only one element of each namelist array, so input
validation rejected residual unused fields before reaching restart. Changing
the fixture to explicit three-element repetition fixed it. Both full checkers
then passed; no solver correction or tolerance change was needed. Existing
compiler warnings remain. No skips in the final checks.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_eb_3d test_reactive_eb_3d_config test_reactive_eb_3d_checkpoint --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_eb_3d_(config|checkpoint)$' --output-on-failure
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 python3 tools/check_eb_ctu_application.py --executable build/spray-amr-reactive-debug/pelef_reactive_eb_3d --work-dir build/eb-inlet-composition-checks/fixed --inflow --inlet-composition
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 python3 tools/check_eb_ctu_application.py --executable build/spray-amr-reactive-selected-debug/pelef_reactive_eb_3d_selected --work-dir build/eb-inlet-composition-checks/selected --inflow --inlet-composition --selected
python3 -m py_compile tools/check_eb_ctu_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, project contract and whitespace passed. Mechanism SHA:
0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4.
Pinned PeleC: bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8.
Pinned PelePhysics: 6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b.

Artifacts: build/eb-inlet-composition-checks/ with build.log,
build-selected.log, LastTest.log, fixed-final.log, selected-final.log.
Final run directories: fixed/run-twa4v2w5 and selected/run-x5uwy3_s.
Initial fixture failures: fixed/run-s27cd_uy and selected/run-hq7cej40,
with fixed.log and selected.log.

## Limits

One overridden inlet face was exercised. Multiple simultaneous mixtures,
other selected mechanisms, single-level MPI/spray, full build matrix,
installed-binary checks, external solver parity and physical validation were
not run. Full PeleC parity remains unproven.
