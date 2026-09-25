# Time-centered external force work — 2026-09-26

Implementation `7567d1f4`, parent `5cb69453`, local main, version 0.246.0.

The full regular-grid split API accepts `external_force(3)` as an alternative
to frozen conserved rates. After the first reaction/transport halves it builds
gas halos and old force/work sources, advances CTU, recomputes work from the
new velocities and adds `0.5*dt*(new_source-old_source)`. Constant force leaves
the momentum increment unchanged. Temperature is recovered after correction.
This matches the old/new source averaging in pinned `Advance.cpp::construct_Snew`
for this contribution, within PeleF's existing split composition; it does not
implement the full reference SDC iteration.

Optional `source_change` returns integrated volume sources, separately from
`boundary_change`, for either frozen or external-force rates. It is published
only after the entire split succeeds. External force requires PPM/PeleC;
simultaneous frozen-source input is rejected. Existing state rollback is retained.

## Verification

`unit_reactive_ppm_step_3d` passed **1/1**, 5.06 s test / 5.11 s total.
Uniform forcing checks exact momentum increment, the quadratic kinetic-energy
term `0.5*dt^2*|force|^2/rho`, unchanged temperature and integrated-source ledger.
A nonuniform outflow case checks final inventory against boundary plus source
changes. Existing source/split/rollback checks also ran. No test failures/skips.
The initial build failed because an automatic array preceded the declaration
of its shape argument; declaration order was fixed. Both build logs are retained.
Provenance, Python syntax, project contract and whitespace passed.
GNU Fortran 13.3.0, Python 3.12.3, existing tests-enabled Ninja Debug build;
existing compiler warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_step_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_step_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/external-force-step-checks/{provenance.log,build.log,build-fixed.log,ctest.log,unit.log}`.
Reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`, `Source/External.cpp` and
`Source/Advance.cpp`, inspected in `/mnt/sn850/Projects/gtt2/external/PeleC`.
Application input, MPI/EB routing, reacting-force qualification and physical
validation remain open.
