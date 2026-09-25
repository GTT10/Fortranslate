# Pinned periodic CTU with flattening and viscosity — 2026-09-26

Comparison implementation `6b3ca71c`, parent `7ecced2b`, local main, version
0.246.0. No numerical solver, tolerance or golden data changed.

The [full periodic CTU oracle](ppm-pinned-ctu-20260926.md) now compiles unchanged
pinned `flatten`, `pc_divu` and `pc_artif_visc` bodies as well as the twelve-solve
pipeline and final species normalization. The array shim supports the reference
bounds checks. Reference conserved states and nodal divergence are formed for
artificial viscosity; the unused temperature flux slot is retained separately.

## Results

**10 cases / 28,800 final-flux values passed**, maximum component-scaled error
**2.988913928686505e-15**, unchanged acceptance bound `3e-12`.
The original six cases remain. Four new cases use a periodic 4x4x4 pressure jump
from 1e5 to 8e5 Pa, temperature 900 K, opposing x velocities +40/-40 m/s and
varying transverse velocities at dt=1e-5. They cover flattening off/on crossed
with artificial viscosity 0/0.1. Both the reference and Fortran output tables
must change when flattening alone or viscosity alone is enabled; inactive
fixtures cannot pass that assertion. Component scaling remains the maximum
absolute face flux per case, floored at one.

The original normal/CTU unit passed **1/1**, 6.41 s test / 6.46 s total.
Provenance, Python syntax, project contract and whitespace passed. No failures
or skips. GNU Fortran/C++ 13.3.0, Python 3.12.3, existing Ninja Debug tests-enabled
`build/spray-amr-reactive-debug`. Build warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_normal_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tools/check_pinned_ppm_ctu.py --reference /mnt/sn850/Projects/gtt2/external/PeleC --executable build/spray-amr-reactive-debug/tests/test_reactive_ppm_normal_3d --work-dir build/ppm-pinned-shock-checks/oracle --shock-viscosity
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_normal_3d$' --output-on-failure
python3 -m py_compile tools/check_pinned_ppm_ctu.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/ppm-pinned-shock-checks/{build.log,comparison.log,ctest.log,unit.log}`
and `oracle/` with pinned sources, compiler log, inputs, output tables, thermo,
source/shim/EOS-adapter/executable hashes and report. Reference:
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

This remains a compiled reference-body comparison with a NASA7 SI adapter,
not a full AMReX/PeleC executable validation. Physical boundaries, WENO, forcing,
EB, selected mechanisms, MPI/application execution, installed audit and physical
shock validation are outside the comparison. Complete PeleC parity is unproven.
