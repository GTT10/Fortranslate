# Selected-mechanism MPI spray application — 2026-09-17

Feature parent `c25d388` (including the separately committed nonnegative chemistry
closure correction), plus this change. Logs and source hashes are retained in
`build/spray-selected-mpi-checks/`.

Both fixed/selected MPI frontends now use a shared mechanism-sized application.
The selected frontend loads the generated thermodynamic/reaction/transport data
through the existing selected validator; checkpoint identity uses the bundle SHA.
Default H2/O2/N2 fractions are assigned by name. The profile-specific selected
MPI 1D frontend is generated only for its supported 2/10-species profiles; other
MPI applications are still generated for larger supported bundles.

A fresh GNU Fortran 13.3.0 Debug / Open MPI 4.1.6 build was configured:

```sh
cmake -S . -B build/spray-selected-mpi-gnu-debug -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran \
  -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec.openmpi \
  -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_SPRAY=ON -DBUILD_TESTING=ON \
  -DMPI_mpi_LIBRARY=/usr/lib/x86_64-linux-gnu/libmpi.so \
  -DMPI_hwloc_LIBRARY=/usr/lib/x86_64-linux-gnu/libhwloc.so \
  -DPELEF_MECHANISM_BUNDLE="$PWD/build/port-methanol-reference/ffcm1_red.json"
cmake --build build/spray-selected-mpi-gnu-debug --target pelef_mpi_spray_3d_selected --parallel 8
```

The 21-species FFCM1_Red bundle SHA is
`ec0b8bfcd71285ba1551c30df0db95c7b23976b71d4184a8025064f7c057f6a8`;
pinned YAML SHA is `ddbba8c97d929530cfb213da45485cca4f979dccf97c1063721bd855d7b3f259`.
`ldd` confirms Open MPI libraries; hwloc resolves to the host oneAPI runtime as
in the preceding MPI environment, so no fully isolated runtime is claimed.

Run twice, with N=2 and N=4:

```sh
python3 tests/spray/test_mpi_spray_application.py \
  --executable build/spray-selected-mpi-gnu-debug/pelef_mpi_spray_3d_selected \
  --mpiexec /usr/bin/mpiexec.openmpi --ranks N --species-count 21 \
  --input cases/mpi_spray_3d/periodic_methanol_water.nml \
  --work-dir build/spray-selected-mpi-gnu-debug/selected-application-N
```

Both selected checker invocations pass, including 21 gas species, strict
nonnegative species, closure, budgets, output refusal, same-rank exact final CSV
comparison, changed-rank continuation and repeated injection at a saved event.
Final selected artifacts are in `selected-application-2/run-z86s01gb` and
`selected-application-4/run-d0qg7j6o` below that build directory. The fixed-mechanism
application CTests also pass 2/2 (15.22 s), no failures/skips, after sharing the
application implementation. Earlier development configure attempts exposed the
old 1D-profile restriction and an unconditional target reference; both were fixed.

Fixed/selected source provenance, project contract, Python syntax and diff checks
pass. No full configuration/installed-executable matrix, all-species-count support
qualification, detailed-fuel physical validation or scaling study is claimed.
