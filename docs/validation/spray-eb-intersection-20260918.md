# Particle segment / nodal EB intersection — 2026-09-18

Implementation: `c80b2165b222cecbbf2faaafcb740dc9686898dc`, based on
`89fac80db0877eedfe03cbe42d4e0689f6ef2417`.

The serial query consumes the original nodal level set and uses the gas EB
builder's six tetrahedra per cell. Barycentric segment clipping and linear
level-set roots find first fluid-to-solid contact, including a solid layer
between fluid endpoints. Normals follow the physical level-set gradient into
fluid, accounting for unequal cell spacing. Invalid inputs preserve all outputs.
Closed-domain endpoints and surface-inward/outward starts have explicit behavior.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja; `/usr/bin/gfortran`,
`/usr/bin/mpifort`, `/usr/bin/mpiexec`. Existing fixed tests-enabled and selected
H2/O2 tests-disabled configurations were reused without configuration changes.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_eb_intersection --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_intersection$' --output-on-failure --output-log build/spray-eb-intersection-checks/ctest-final.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_spray_tree_runtime --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_mechanism_provenance.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **1/1 passed**, test time **0.00 s**, total **0.04 s**. The test checks an
analytic oblique plane on anisotropic cells, a fluid-to-fluid segment crossing
a solid slab, surface starts, a solid-start rollback, closed upper endpoints,
a stationary fluid point and exterior rejection. Both library builds passed;
the selected configuration was compile-checked only. Provenance, Python syntax,
project contract and whitespace checks passed. No build/test failures or skips.
Compiler real-comparison warnings remain. An earlier successful focused run
preceded the final integer-dimension overflow guard.

Artifacts: `build/spray-eb-intersection-checks/` contains `build-final.log`,
`ctest-final.log`, `build-selected.log`, `build-selected-final.log` and initial
build/test logs. The first target in the selected build directory was the fixed
library; the second explicitly checked the selected library.

This is a bounded geometry primitive (at most one million cells; swept cell-box
scan), not an integrated EB particle solver. It does not supply EB gas sampling,
deposition, impact response, film evolution, MPI ownership or AMR stepping.
No external solver parity, convergence, production scaling or physical
validation was performed.
