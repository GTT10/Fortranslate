# MPI AMR reacting-spray step — 2026-09-18

Implementation: `f0ad69128dd3bb278795703c1580468eb54899b5`, based on
`e7d2eae`. `advance_mpi_amr_spray_reactive` composes particle half-intervals
around the existing reaction/transport/hydro split. Gas updates may use common
or recursive substeps. The caller selects the overall timestep and retry policy;
particle halves retain their collective adaptive substeps. U/T, parcels,
component masses, optional ID watermark and gas correction/limiter diagnostics
commit only after all stages succeed. Failure returns zero particle substeps
and preserves caller state and diagnostics.

The MPI spray tree runtime now links the existing MPI gas runtime, removing
its duplicate cell-query/restriction/thermodynamic source compilation. The
selected tree runtime links the selected MPI gas runtime, preserving the
mechanism boundary. Periodic faces are used throughout. Executable setup,
injection, regrid and checkpoint orchestration remain open.

## Completed focused verification

Clean GNU Fortran 13.3.0 Debug/Ninja builds used `/usr/bin/gfortran`, Open MPI
`/usr/bin/mpifort` and `/usr/bin/mpiexec`. Fixed tests were enabled; selected
H2/O2 tests were disabled. No installed prefix was used.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/spray-amr-reactive-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_SPRAY=ON -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=ON
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_reactive --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_reactive_(2|4)$' --output-on-failure --output-log build/spray-amr-reactive-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/spray-amr-reactive-selected-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_SPRAY=ON -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=OFF -DPELEF_MECHANISM_BUNDLE=mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, no failures/skips, **75.45 s** total: two ranks **38.16 s**,
four ranks **37.25 s**. The three-level water/peroxide fixture starts all
particles on one rank. An invalid gas Riemann solver rejects after the first
particle half, preserving gas, temperature, particle position/velocity/mass/
temperature, component inventory, ID watermark and diagnostic sentinels.

Both common and recursive gas paths complete with viscosity, conduction and
species diffusion enabled, PCM reconstruction and HLLC fluxes. Composite gas
plus liquid mass, momentum, energy and species inventory differences remain
within 3e-12 of the respective inventory scale. A subsequent implicit-chemistry
step succeeds and preserves total mass, momentum and energy within the same
bound. This is a coupling/conservation fixture, not ignition validation or
an independent accuracy/order oracle for the split scheme.

The selected library compiled. Its Ninja dependency edge includes the selected
MPI tree runtime and excludes `pelef_core`; this is compile/dependency evidence,
not selected executable qualification. Provenance, Python syntax, project
contract and whitespace checks passed. No build/test attempt failed; compiler
warnings remain. No Python source changed. No full matrix, install audit,
external solver parity, AMR convergence, scaling or physical spray validation
was performed.

Logs: `build/spray-amr-reactive-checks/{configure,build,ctest,configure-selected,build-selected}.log`.
Executable: `build/spray-amr-reactive-debug/tests/test_spray_mpi_amr_reactive`.
