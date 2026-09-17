# Distributed AMR particle field deposition — 2026-09-18

Implementation: `a594cd156a8454ff551413b60fb526f264eb6eb8`, based on
`9332646`. `deposit_mpi_amr_tree_spray_field` shares composite stencil preparation
with MPI sampling, packs destination-owner cell coordinates and volume-scaled
contributions, and exchanges them directly with Alltoallv. Each owner applies
source ranks in increasing order while retaining source particle/stencil order.
The existing MPI bottom-up restriction refreshes covered parent cells. Only a
fully accepted private candidate replaces the original owner-local state.

Outgoing request and incoming contribution capacities are independently bounded
per rank (each defaults to one million). Nonfinite owned fields/contributions,
metadata/stencil failures, capacity limits and unsafe density/addition arithmetic
reject collectively. Globally empty feedback validates inputs and preserves the
field, including any prior covered-parent values. No full-field gather is used.
Changed particle ownership may change accumulation roundoff. This generic field
API does not yet enforce EOS admissibility or recover MPI temperatures.

## Completed focused verification

Clean GNU Fortran 13.3.0 Debug/Ninja builds, `/usr/bin/gfortran`, Open MPI
`/usr/bin/mpifort` and `/usr/bin/mpiexec`: fixed tests enabled and selected H2/O2
library tests disabled. No installed prefix is involved.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/spray-amr-deposition-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_SPRAY=ON -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=ON
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-deposition-debug --target test_spray_mpi_amr_deposition test_spray_mpi_amr_sampling --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-deposition-debug -R '^(spray_mpi_amr_deposition_(2|4)|spray_mpi_amr_sampling_2)$' --output-on-failure --output-log build/spray-amr-deposition-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/spray-amr-deposition-selected-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_SPRAY=ON -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=OFF -DPELEF_MECHANISM_BUNDLE=mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-deposition-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**3/3 passed**, no failures/skips, **1.45 s** total: deposition on two ranks
**0.46 s**, existing sampling on two ranks **0.47 s**, deposition on four ranks
**0.48 s**.

Three-level, two-component fixtures compare every owned cell, including averaged
parents, exactly with serial deposition in matching source order. They exercise
coarse/fine, nested and periodic clouds, multiple source ranks contributing to
the same cells, a zero-particle gas owner, a rank with no gas storage, changed
gas owner maps and globally empty feedback. A separate exposed-cell volume
sum followed by MPI reduction matches the initial composite inventory plus
signed particle contributions within 3e-13 of the inventory scale. Receive
capacity, nonfinite feedback and density-overflow failures preserve every owned
cell. Nonfinite initial fields reject. Sampling regression covers the shared
stencil preparation refactor and its metadata/position failures.

The selected library compiled successfully; its Ninja dependency edge contains
only selected runtimes, without `pelef_core`. This is compile/dependency evidence,
not selected executable or installed-binary evidence. Provenance, Python syntax,
project contract and whitespace checks passed. No build/test attempt failed;
compiler warnings remain. No Python source changed. No allocation-failure
injection, MPI fault recovery, scaling, gas EOS coupling, AMR convergence or
physical spray validation is claimed.

Logs: `build/spray-amr-deposition-checks/{configure,build,ctest,configure-selected,build-selected}.log`.
Executables: `build/spray-amr-deposition-debug/tests/{test_spray_mpi_amr_deposition,test_spray_mpi_amr_sampling}`.
