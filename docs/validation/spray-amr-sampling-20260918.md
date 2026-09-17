# Owner-local MPI AMR particle sampling — 2026-09-18

Implementation: `22093f816198a8df381d70d85a552892fd4de36c`, based on
`fa4c3d3`. `sample_mpi_amr_tree_spray_field` builds each rank's local composite
cloud stencils and uses the existing sparse tree-cell fetch API to request
only their exposed cells. Values accumulate in serial stencil order. The full
gas field is never gathered. The generic cell-query module is compiled into
the isolated fixed/selected spray-tree runtime without a fixed-mechanism
dependency in the selected target.

Topology, owners, component count, physical bounds/spacing, periodic flags and
local field allocation/shapes are checked collectively. Only requested payload
values require finiteness; covered/unrequested values are not scanned. Invalid
positions, metadata, requested values, allocation or arithmetic rejection
preserve all ranks' prior sample arrays. Each rank may bound total requests
independently (zero to one million, default one million). Zero-query ranks still
serve requests from other ranks. This adds sampling, not MPI deposition/EOS
feedback or a coupled particle time integrator.

## Completed focused verification

Clean GNU Fortran 13.3.0 Debug/Ninja builds, `/usr/bin/gfortran`, Open MPI
`/usr/bin/mpifort` and `/usr/bin/mpiexec`: fixed tests enabled and selected
H2/O2 library tests disabled, without any installation prefix.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/spray-amr-sampling-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_SPRAY=ON -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=ON
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-sampling-debug --target test_spray_mpi_amr_sampling --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-sampling-debug -R '^spray_mpi_amr_sampling_(2|4)$' --output-on-failure --output-log build/spray-amr-sampling-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/spray-amr-sampling-selected-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_SPRAY=ON -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=OFF -DPELEF_MECHANISM_BUNDLE=mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-sampling-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, no failures/skips, **0.99 s** total: two ranks **0.46 s**,
four ranks **0.49 s**. Three-level, two-component fields match serial samples
exactly at coarse/fine, nested and periodic boundaries. The fixture includes
unequal local query counts, a zero-query gas owner, a rank with no gas ownership,
changed owner maps and globally empty requests. Unrequested covered-cell NaNs
do not affect samples; requested NaNs, inconsistent bounds/owners, a missing
owned field, a bad position on one rank and one rank's request-capacity failure
all reject collectively without changing samples.

The selected library compiled successfully; its Ninja dependency edge uses only
selected runtime libraries, without `pelef_core`. This is compile/dependency
evidence, not a selected executable or installed-binary audit. Provenance,
Python syntax, project contract and whitespace checks passed. No build/test
attempt failed; compiler warnings remain. No Python source changed. No memory
fault injection, MPI fault recovery, scaling, AMR convergence or physical spray
validation is claimed.

Logs: `build/spray-amr-sampling-checks/{configure,build,ctest,configure-selected,build-selected}.log`.
Executable: `build/spray-amr-sampling-debug/tests/test_spray_mpi_amr_sampling`.
