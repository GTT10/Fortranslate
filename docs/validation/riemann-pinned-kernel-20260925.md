# Compiled pinned primitive Riemann kernel — 2026-09-25

Source: base `4421987`, comparison tool `9cdbac0`. GNU Fortran/C++ 13.3.0,
Ninja Debug, Python 3.12.3; reused fixed tests-enabled
`build/spray-amr-reactive-debug`. Artifacts: `build/riemann-pinned-kernel-checks/`.

The comparison tool compiles the unchanged `riemann` function extracted with
`git show` from PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
Riemann.H SHA256:
`073ae14b2f7925de2858ef312d02fc50f033c2a344dccc5b24ecde1f77ed8253`.
The C++ driver supplies SI NASA7 energy/sound-speed methods using molecular
weights and low/high coefficients exported by the Fortran test executable.
This independently evaluates EOS arithmetic with shared thermodynamic data;
it does not compile or validate the complete PelePhysics EOS implementation.
The adapter supplies pinned stationary and wave-regularization constants.
Rotating-frame inputs are zero and boundary velocity multiplier is one.

The Fortran stream calls the production primitive-state Riemann entry directly.
Outputs are all 12 conserved flux components, interface normal velocity and
interface pressure. The reference retains its species flux products; PeleF
retains final-species closure by subtraction. These arithmetic differences are
included in the comparison, not masked by rewriting reference output.

## Completed checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target test_reactive_pelec_riemann --parallel 2`: passed.
- With those library variables unset:
  `python3 tools/check_pinned_riemann_kernel.py --reference /mnt/sn850/Projects/gtt2/external/PeleC --executable build/spray-amr-reactive-debug/tests/test_reactive_pelec_riemann --work-dir build/riemann-pinned-kernel-checks/oracle`: passed.
  **56 cases / 784 values**; maximum `abs(a-b)/max(1,abs(a),abs(b))`
  **3.111830054382345e-14**, fixed limit `3e-12`.
  Cases cover equal states, material contacts and unequal pressure/density/
  composition, velocities `-2500`, `-100`, `-1e-8`, `0`, `1e-8`, `100`, `2500`
  m/s with bounded perturbations, temperatures 600–1800 K, pressure around
  100 kPa and nonzero transverse velocity. Random seed `20260925`.
  An initial 49-case run also passed before adding explicit material contacts.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset,
  `HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`:
  `ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_pelec_riemann$' --output-on-failure`:
  **1/1 passed**, 0.00 s test / 0.05 s total; no failures/skips.
- `python3 -m py_compile tools/check_pinned_riemann_kernel.py tools/check_pinned_ppm_kernels.py tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, `git diff --check`: passed.

The artifact `oracle/` contains original reference header, extracted source,
exported thermodynamic coefficients, common inputs, both output tables and a
report with executable hash and exact C++ compile command. No solver change was
required by this comparison. Full PeleC executable/EOS parity, normal projection,
complete CTU, convergence, full configuration matrix and physical validation
were not run. This is bounded kernel evidence only.
