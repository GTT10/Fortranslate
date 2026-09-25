# PPM CTU source routing — 2026-09-26

Implementation `21e0bdfa`, parent `29a88e75`, local main, version 0.246.0.

Periodic and supplied-halo PPM CTU flux APIs now accept paired primitive and
internal-energy-density source arrays. The periodic wrapper fills four source
ghost layers; the halo API requires caller-supplied matching layers. Sources
require the PeleC pressure policy. Both final double-transverse face states
receive half-step sources, after conservative conversion, matching `pc_transdd`.
Single-transverse corner solves and final inflow PLM replacement retain their
reference ordering. Missing pairs, wrong shapes and nonfinite rates fail before
publishing fluxes. Absent optional arrays retain the prior path.

## Verification

The full compiled-reference oracle passed **28 cases / 80,640 final-flux values**,
maximum component-scaled error **2.988913928686505e-15**, unchanged bound `3e-12`.
The existing eighteen cases remain. Ten added cases use independent mixtures,
flattening enabled, dt=1e-5 and artificial viscosity 0/0.1, crossing periodic,
outflow, slip, no-slip and inflow boundaries with finite constant sources.
Source rates are rho=0.2, velocity=(300,-400,200), pressure=-20000,
first two species=(10,-10), and internal-energy density=5000, in corresponding
SI units per second. Other species sources are zero.

Related normal and halo units passed **2/2**, 6.33 and 6.54 s, 12.92 s total.
No failures or skips. Provenance, Python syntax, contract and whitespace checks
passed. GNU Fortran/C++ 13.3.0, Python 3.12.3, existing tests-enabled Ninja Debug
build; existing compiler warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_normal_3d test_reactive_ppm_halo_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tools/check_pinned_ppm_ctu.py --reference /mnt/sn850/Projects/gtt2/external/PeleC --executable build/spray-amr-reactive-debug/tests/test_reactive_ppm_normal_3d --work-dir build/ppm-ctu-source-checks/oracle --shock-viscosity --physical-boundaries --sources
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(unit_reactive_ppm_normal_3d|unit_reactive_ppm_halo_3d)$' --output-on-failure
python3 -m py_compile tools/check_pinned_ppm_ctu.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/ppm-ctu-source-checks/` holds provenance, build, comparison,
CTest and unit logs; `oracle/` retains pinned sources, generated driver, inputs,
outputs, compiler log and report with source/shim/adapter/executable hashes.
Reference commit: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

This compares extracted reference bodies with the NASA7 SI adapter. Spatially
varying sources, conserved-source conversion, application source generation and
time integration, MPI source exchange, EB, selected mechanisms and physical
validation are outside this check. Full PeleC equivalence remains unproven.
