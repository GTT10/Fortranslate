# Shock-flattening denominator correction — 2026-09-25

Code `1a36c4d`, parent `b3688ea`, changes both pressure-difference denominator
floors in `reactive_ppm_flattening_coefficient` from machine `tiny` to `1e-200`.
This matches `small_pres` in pinned PeleC
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`, `Source/Godunov.H:208-259`.
The sensor is shared by regular PPM/PLM and EB CTU preparation. With a zero
outer pressure difference and an ordinary nonzero inner jump, dividing by
machine tiny raised overflow before coefficient saturation. EOS temperature
limits, shock thresholds and output clipping are unchanged.

## Focused evidence

GNU Fortran 13.3.0 / OpenMPI 4.1.6 / Python 3.12.3, existing Ninja Debug build.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_halo_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_halo_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The added regression failed against the old implementation with
`IEEE_OVERFLOW_FLAG` (0.11 s test / 0.15 s total). After the correction,
**1/1 passed**, 6.40 s / 6.44 s total. Eight sensor cases cover direct/shifted
zero denominators, positive/negative pressure jumps, and compression/expansion.
They check both the expected coefficient (zero or one) and absence of the
IEEE overflow flag. Existing regular halo/CTU checks also pass in that program.
No additional broad matrix was run; this does not assert overflow immunity for
arbitrarily large pressure inputs or establish external CFD parity.

Artifacts: `build/flattening-denominator-checks/{before-build.log,after-build.log,LastTest-before.log,LastTest-after.log}`.
No build failures or skips. The pre-fix regression failure is intentional
reproduction evidence. Compiler warnings remain. Provenance matched
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`;
Python syntax, project contract and diff checks passed. No golden data or
verification tolerances changed.
