# Compiled pinned CTU transverse kernels — 2026-09-26

Comparison implementation `b5a60e39`, parent `7dd04b28`, local main, version
0.246.0. No solver arithmetic, tolerances or golden data changed.

`tools/check_pinned_ppm_transverse.py` reads Godunov.H at manifest-pinned PeleC
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8` and compiles unchanged bodies of
`pc_transdo_passive`, `pc_transdd_passive`, `pc_transdo` and `pc_transdd`.
The small C++ shim supplies array indexing, integer geometry and numeric
constants. No EOS adapter is needed. Fortran stream mode calls the production
`reactive_ppm_transverse_pressure_3d` helper with the same conserved face state,
weighted flux differences, Godunov pressure/velocity and gamma.

## Results

**144 cases / 2,304 values passed**, maximum scaled error
**1.4551886041586193e-16**, fixed limit `3e-12`. Inputs combine one/two transverse
directions, three normal axes, regular/signed-species/pressure-floor modes and
eight random samples (seed 20260926). Both output sides are compared for density,
three velocities, pressure, two species weights and internal energy. Both sides
receive the same initial state. Pressure-floor outputs additionally require exact
agreement, avoiding a vacuous absolute-tolerance check on tiny values.

The unchanged ordinary halo unit passed **1/1**, 6.40 s test / 6.45 s total.
Provenance, Python syntax, project contract and whitespace checks passed.
No failures or skips occurred. GNU Fortran/C++ 13.3.0, Python 3.12.3, existing
Ninja Debug tests-enabled `build/spray-amr-reactive-debug`. Build warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_halo_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tools/check_pinned_ppm_transverse.py --reference /mnt/sn850/Projects/gtt2/external/PeleC --executable build/spray-amr-reactive-debug/tests/test_reactive_ppm_halo_3d --work-dir build/ppm-pinned-transverse-checks/oracle
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_halo_3d$' --output-on-failure
python3 -m py_compile tools/check_pinned_ppm_transverse.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts under `build/ppm-pinned-transverse-checks/`: `build.log`,
`comparison.log`, `ctest.log`, `unit.log`, and `oracle/` with source, generated
C++ driver, input/output tables, compiler command and source/executable hashes.

These are regular-cell, zero-source kernel checks with two passive species.
They do not test EB aperture gates, forcing, additional passive groups, EOS,
normal prediction or the ordering/assembly of the twelve-solve CTU pipeline.
Full-CTU external parity, application/MPI runs, install audit, full matrix and
physical validation remain outside this evidence.
