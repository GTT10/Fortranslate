# Regular spray PLM coefficient input/restart — 2026-09-25

## Source and scope

Code `0c51be1`, parent `c4f2336`. Serial spray accepts explicit nonnegative
PCM/PLM coefficients; regular MPI spray accepts PLM coefficients (its supported
schemes remain PLM/PPM). Both already forward the coefficient into hydro.
Automatic -1 remains 0.1 for PPM/pelec and zero otherwise.

Positive non-PPM coefficients are now recorded in application checkpoint
context. Zero/automatic non-PPM retains the legacy context. Existing PPM
serialization is unchanged. Changed coefficients reject restart. MPI retains
collective coefficient agreement. The checker no longer expects non-PPM
coefficient input to reject before restart; reconstruction changes are checked
by checkpoint identity instead.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Ninja Debug, Python 3.12.3. Reused fixed
tests-enabled `build/spray-amr-reactive-debug`. Focused CTest **2/2 passed**:
serial 28.18 s, MPI 8.85 s; total 37.07 s. No failures or skips. Compiler
warnings remain.

Serial PLM uses coefficient 0.05 with Mach inflow through coupled stepping.
Restart matches uninterrupted output; zero coefficient changes gas output.
Changed coefficients, -2 and NaN reject without final history output.
Automatic -1 and explicit zero continuation CSVs match exactly.

MPI PLM uses coefficient 0.037 in the periodic mixture fixture, including
same-rank and 2-to-4-rank continuation. Existing state/particle/history
comparisons pass. Changed coefficients and invalid/nonfinite inputs reject;
automatic/zero changed-rank continuation records match. The coefficient fixture
uses the first accepted reference time as its split point, as already done by
the PPM coefficient test.

This is fixed-mechanism coupled application/restart evidence. Selected builds,
serial PCM runtime, non-PPM AMR spray restart, a full matrix and physical
validation were not independently run here. Mechanism provenance, Python
syntax, project contract and whitespace checks passed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_spray_3d pelef_mpi_spray_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_(mpi_)?plm_difmag_restart$' --output-on-failure
python3 -m py_compile tools/check_spray_boundaries.py tests/spray/test_mpi_spray_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/spray-plm-difmag-input-checks/build.log`, copied
`LastTest.log` and fixture paths recorded there.
