# MPI EB ancestor time interpolation — 2026-09-18

Implementation: `4d521eae5dc324f7e62d87dba16da069f69fb30e`, based on
`448f4a6`. Extends `fill_mpi_eb_tree_halo_3d` with optional owner-local
`ancestor_start`, `ancestor_end`, `temperature_start`, `temperature_end` and
`ancestor_fraction` arguments.

All five arguments must be present together on every rank. Endpoint containers
use the supplied topology, ownership and immutable geometry. Fractions have
length equal to the maximum tree level, indexed root first, must be finite in
[0,1], and must agree across ranks. Every donor below the target level uses its
own depth's interpolated conserved state and temperature guess, followed by
EOS recovery. Endpoint fractions select the corresponding endpoint exactly.
Same-level peers and physical outflow retain the current averaged snapshot.

Provisional endpoint fields are not averaged from descendants. Interpolation
uses owner-local copies; only requested ancestor cells are exchanged. No full
state gather is introduced. Output U/T halos publish together after collective
success. Inputs remain unchanged. This supplies the temporal donor operation;
MPI recursive stepping is still separate work.

## Completed focused verification

Fresh GNU Fortran 13.3.0 Debug/Ninja build, compiler `/usr/bin/gfortran`,
Open MPI wrapper `/usr/bin/mpifort`, launcher `/usr/bin/mpiexec`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/mpi-eb-tree-temporal-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=ON
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/mpi-eb-tree-temporal-debug --target test_mpi_eb_tree_halo_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/mpi-eb-tree-temporal-debug -R '^unit_mpi_eb_tree_halo_3d_[23]$' --output-on-failure --output-log build/mpi-eb-tree-temporal-checks/ctest.log
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, no failures/skips, total **1.73 s** (two ranks 0.88 s,
three ranks 0.81 s). A five-node, three-level cut-cell tree exercises both
distributed and single-owner layouts, including empty ranks. Fractions [0,1],
[1,0] and [0.25,0.75] produce exact serial U/T halos including inactive zeros.
Endpoint states differ by depth and remain unmodified. Existing spatial halo,
geometry mismatch and current-donor EOS checks also pass.

Added rejection cases cover rank-disagreeing fractions, partial context,
rank-disagreeing context presence and invalid interpolated ancestor density.
Every rejected call preserves previous halo outputs. An unused endpoint filled
with NaNs is ignored at fraction zero.

Provenance, configure/build, Python syntax, project contract and whitespace
checks pass. Compiler warnings remain; no failed build/test attempt occurred.
No Python source changed. This is a focused API check, not an application
subcycling, temporal convergence, installed-binary, scaling, full-matrix or
physical-validation result.

Logs: `build/mpi-eb-tree-temporal-checks/{configure,build,ctest}.log`.
