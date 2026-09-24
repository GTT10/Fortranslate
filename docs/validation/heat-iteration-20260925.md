# Thermal Spalding iteration — 2026-09-25

Implementation `deea9ab`, base `52d2523`. Reference PelePhysics
`6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`, `Drag.H:9–46`,
SHA256 `bcc2bcb2f01b8cedfeb92d0bfc9ae2c8e723e77285f7011aeb77d5b138a29345`,
verified against `references/spray_property_scope.json`.

The AS correction now initializes B_T using Nu_0, directly updates B_T through
the corrected Nu, stops at absolute B_T change <=1e-7, and recomputes final Nu
from the final B_T. This replaces damped Nu iteration. The initial correction
and maximum 100 subsequent updates follow pinned `calcHeatCoeff`.
PeleF retains failure on nonconvergence, an exponent guard, its B_M<=1e-12
uncorrected limit and small-B blowing-factor limit. The reference instead
returns its last iterate and has a B_M<=1e-15 zero-return branch.
Binary diffusivity and temporal integration remain different conventions.
This establishes focused positive-B coefficient agreement, not full solver parity.

Numerical restart identifiers changed to `spray-v4-heat-iteration` and
`mixed-liquid-components-v6-heat-iteration` so old model contexts are rejected.

## Independent reference and checks

The pinned C++ function was extracted unchanged into a standalone harness
(only an `amrex::Real=double` alias and standard headers supplied). At B_M=0.25,
Pr=0.6, Sc=0.72, fuel Cp=2000 J/(kg K), rho=1 kg/m3, D=2.5e-5 m2/s,
lambda=0.03 W/(m K), mu=1.8e-5 Pa s and diameter=1e-4 m, it returns:

| Slip m/s | Corrected heat Nu |
| --- | --- |
| 0 | 1.6510934499401388 |
| 1 | 2.2947576817546693 |
| 50 | 7.2882103277227754 |
| 500 | 15.03362047582373 |

`spray_abramzon` now compares all four against the Fortran coefficient with
relative tolerance 2e-13. Existing tolerances were unchanged.

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug, Python 3.12.3; existing fixed
tests-enabled and selected H2/O2 tests-disabled builds. Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
c++ -O0 build/heat-iteration-checks/oracle.cpp -o build/heat-iteration-checks/oracle
build/heat-iteration-checks/oracle
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_core pelef_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_(abramzon|vapor_cp|mixture_rates|heat)$' --output-on-failure --output-log build/heat-iteration-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/spray/test_spray_application.py --exe build/spray-amr-reactive-debug/pelef_spray_3d SprayApplication.test_molecular_weights_restart
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/spray/test_spray_application.py --exe build/spray-amr-reactive-selected-debug/pelef_spray_3d_selected SprayApplication.test_molecular_weights_restart
python3 -m py_compile tests/spray/test_spray_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Core **4/4 passed**, **0.05 s**. Fixed/selected application methods each passed
**1/1**, **1.672/1.662 s**, covering scalar/port/mixture molecular weights,
conservation, exact CSV restart and changed/invalid input rejection on a
2-by-2-by-2 chemistry-disabled case. Builds, provenance, Python syntax,
project contract and whitespace passed. No failures or skips; compiler warnings
remain. Artifacts: `build/heat-iteration-checks/` holds oracle.cpp/oracle.txt,
build-fixed.log/build-selected.log, ctest.log/LastTest.log and
application-fixed.log/application-selected.log.

No MPI execution, full configuration matrix, installed audit, external solver
run or physical validation was performed in this increment.
