# Recursive periodic EB films — 2026-09-27

Implementation `72051c20`, based on `2529e2c3b20fe732d133b14a0dc0c1bbd98750c7`,
local main, version 0.246.0.

## Change and scope

Optional periodic axes now reach film physical-level initialization, virtual and
ghost construction, generating-level newborn attachment and final collection.
Migration collectively checks the axes and matching opposite nodal geometry.
Ownership lookup wraps coordinates; film payload coordinates remain unchanged.
The recursive thermal/source path uses these axes, and application input now
accepts periodic wall impacts with recursive SDC.

Surface-validation probes must still remain inside the closed physical domain.
Normals whose probes cross a periodic seam remain unsupported. Copy selection
retains one destination per input, not multiple images in one population.
This increment does not establish arbitrary periodic EB geometry, physical
validation or exact PeleC parity.

## Independently completed checks

Artifacts: `build/spray-recursive-periodic-film-checks/`. Existing fixed-mechanism
Ninja Debug build, GNU Fortran 13.3.0, OpenMPI 4.1.6, `/usr/bin/mpifort` and
`/usr/bin/mpiexec`. H2O2 source SHA:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc_recursive_(impact|impact-periodic)_2$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc_recursive_periodic_impact_app$' --output-on-failure
python3 -m py_compile tests/application/check_recursive_eb_spray_impacts.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final selected tests: **3/3 passed**, no skips. Nonperiodic and periodic library
impact checks took 8.61 and 8.70 s (17.37 s total); the periodic application check
took 36.74 s (36.79 s total). Provenance, builds, Python syntax, contract and
whitespace checks passed. Existing compiler warnings remain.

The library fixture includes mixed pre-existing films, splash/newborn films,
chemistry, transport, thermal exchange and rejection/retry rollback. Its periodic
combined inventory residual was **2.8115820151083464e-15**, with four particles
and three films. The nonperiodic residual was **1.5972576137871884e-15**.

The application crosses the z seam before hitting the planar EB and checks that
the deposited film lies near the lower z face. Two and four ranks each finish
with four particles and one film; both budget residuals are
**1.4634730538842126e-15**. Gas, particle and film outputs match within the
unchanged `3e-12` comparison criterion. Cases:
`build/spray-amr-reactive-debug/spray-sdc-recursive-periodic-impact-app/impact_gc46li9z/`.

Initial development failures: a Fortran line exceeded 132 columns (wrapped);
the first application CMake target name was incorrect (corrected). The first
application fixture made seam arrival coincide with EB impact and deposited at
the valid upper endpoint, failing the check for a lower-face deposit. Moving
initial z from 0.003999999 to 0.0039999995 makes seam arrival precede impact;
the solver and tolerances were not changed for this failure. The failed test
log is retained as `app-ctest-initial.log`.

No periodic film restart, selected-mechanism configuration, full matrix or
external-solver comparison was run in this increment.
