# MPI EB particle feedback and EOS recovery — 2026-09-18

Implementation: `3eed8a4537969cdf1dedfd144ca01516e553a25e`, based on `ce9ebd9`.

Field and gas deposition accept paired optional EB geometry/nodal arguments.
The existing Alltoallv contribution routing uses EB stencil weights and physical
fluid volumes. Bottom-up parent restriction averages the weighted conserved
field and fluid fraction together, then divides only in hidden parent cells.
Restricted fluid fractions must match parent metrics within 128 epsilon relative
to the fraction. Exposed cells retain the exact deposited values; solid cells
retain incoming values. The gas wrapper validates/reconstructs temperatures only
in cells with positive fluid volume and publishes state and temperature together
after collective success. Regular AMR remains the default path.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja; existing fixed tests-enabled
and selected H2/O2 tests-disabled configurations reused unchanged.
Compiler `/usr/bin/gfortran`, MPI wrapper `/usr/bin/mpifort`, launcher
`/usr/bin/mpiexec`.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_deposit --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_deposit_[24]$' --output-on-failure --output-log build/spray-mpi-eb-deposit-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_gas --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_gas_2$' --output-on-failure --output-log build/spray-mpi-eb-deposit-checks/ctest-regression.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tools/check_project_contract.py tools/check_mechanism_provenance.py
python3 tools/check_project_contract.py --source .
git diff --check
```

New MPI checks **2/2 passed**, two ranks **0.48 s**, four ranks **0.51 s**, total
**1.04 s**. Regular AMR gas regression **1/1 passed**, **0.53 s**, total **0.56 s**.
Selected compilation, provenance, syntax, contract and whitespace checks passed.
No failures or skips; existing compiler warnings remain. Selected execution
and the full configuration matrix were not run.

The 4-cubed root and ratio-two child have different owners and a y=1.25 planar
cut. Initial active cells hold N2 at 1000 K and nonzero momentum. Solid cells
hold zero conserved state and sentinel temperature 42. A query at (1.9,1.3,1.5)
deposits across the cut coarse/fine interface. All composite components match
the analytic fluid volume 44 plus a 0.01-volume-equivalent contribution at
relative 1e-13 (scale floor 1). Hidden-parent integrals match fine integrals at
the same tolerance. Solid values remain exact, proving they bypass EOS recovery.
A subsequent large negative contribution fails EOS recovery and leaves every
owner-local gas/temperature array unchanged on all ranks.

Artifacts: `build/spray-mpi-eb-deposit-checks/` contains `build.log`,
`build-selected.log`, `build-regression.log`, `ctest.log`, `ctest-regression.log`.

This connects distributed field feedback, not particle time integration or a
distributed EB spray application. No small-cell source redistribution is added;
the caller must reduce failed source intervals. Geometry remains replicated.
No convergence, production scaling, external parity or physical validation ran.
