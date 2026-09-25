# EB periodic peer flux correction — 2026-09-26

Source `7055ecc`, parent `cc7544d`, local main. The existing transactional
`reconcile_reactive_eb_peer_flux_3d` API adds optional `periodic_length(3)`.
Exactly one component must be positive, with zero tangential components.
The supplied integer patch offset describes adjacency after translation by
one domain length. The physical origin check accounts for that translation;
spacing, aperture and periodic face-centroid compatibility remain required.
The receiver's time-integrated flux and corrected state/temperature publish
together only after EOS recovery succeeds. Ordinary peer behavior is unchanged.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug, Python 3.12.3; existing tests-enabled build.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_reactive_reflux_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_reactive_reflux_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final result **1/1 passed**, 0.01 s test / 0.06 s total, no skips. The peer
cases now run for both ordinary and translated periodic interfaces, all
three directions, both sides and partial tangential overlap. An independently
summed open-face transfer agrees with the volume-weighted inventory change
within the existing 3e-12 tolerance. Matching flux correction is idempotent;
nonphysical correction and negative period preserve state, temperature and
flux arrays. Existing coarse/fine restriction and reflux checks also pass.

Initial test setup failed while constructing translated receiver states.
Moving translation after initialization exposed a second fixture issue:
absolute cut-boundary centroids must move with the geometry, while regular
and covered cells retain zero boundary data. Correcting those fixture
coordinates produced the passing run without relaxing production validation.
Provenance, syntax, project contract and whitespace checks passed.

Artifacts: `build/eb-periodic-peer-flux-checks/{build.log,LastTest.log,
initial-failure.log,geometry-failure.log}`. Mechanism SHA:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

This qualifies the correction kernel only. Serial/MPI tree periodic traversal,
coarse/fine seam synchronization and periodic time updates remain pending.
No full matrix, installed audit, external solver parity or physical validation
was run.
