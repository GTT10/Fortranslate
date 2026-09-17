# MPI EB tree automatic common intervals — 2026-09-18

Implementation: `c039f1d40651a91de50485cb1dee7a1e21dbd2e6`, based on
`c47bf21`. `mpi_eb_tree_timestep_3d_mod` belongs to `pelef_mpi_support` and
is called by the shared application when MPI automatic stepping is enabled.

The estimator averages a private owner-local snapshot, scans owned fluid
patches for wave-speed and diffusivity maxima, then shares those maxima to
cover remote peer/ancestor donors. Each owner bounds its patch Euler and
transport updates, including hidden parents; a collective minimum gives the
common interval. Split transport retains the conservative full-interval bound.
Empty ranks participate without allocating remote field payloads. CFL values
and transport switches must agree across ranks. Rejected estimates preserve
both input containers and the caller's timestep. All-solid trees are rejected.

The application now accepts fixed or automatic common intervals on static
MPI trees, retaining the time-step cap, retries and checkpoint policy fields.
Subcycling and regridding/adaptation remain unsupported in this entry point.
No new checkpoint schema or root field gather is used for interval selection.

## Completed focused verification

Fresh GNU Fortran 13.3.0 Debug/Ninja build; compiler `/usr/bin/gfortran`, MPI
wrapper `/usr/bin/mpifort`, launcher `/usr/bin/mpiexec`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/mpi-eb-tree-auto-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=ON
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/mpi-eb-tree-auto-debug --target test_mpi_eb_tree_transport_3d pelef_reactive_eb_tree_3d pelef_mpi_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/mpi-eb-tree-auto-debug -R '^(unit_mpi_eb_tree_transport_3d_[23]|regression_mpi_eb_tree_auto_3d)$' --output-on-failure --output-log build/mpi-eb-tree-auto-checks/ctest.log
python3 -m py_compile tools/check_mpi_eb_tree_application_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**3/3 tests passed**, zero failures/skips, total **26.75 s**. Two/three-rank
unit tests took 4.02/3.91 s; the application regression took 18.79 s.

The existing transport fixture now checks hydro-only, coupled hydro/diffusion
and diffusion-limited estimates against the serial estimator. Three cut-plane
orientations, four/five patches, distributed/single-owner layouts and two/three
ranks give **72 exact estimate comparisons**, including empty ranks. Estimate
calls leave U/T unchanged; rank-inconsistent CFL values preserve sentinel dt.
Existing transport state/conservation/rollback checks also pass.

The application uses `hydro_cfl=5e-5`, starts with
**6.8172036428352156e-10 s** below its 1e-9 s cap, and reaches 2e-9 s in three
steps. The final clipping gives minimum interval 6.3655930903890169e-10 s.
Serial, two-rank, three-rank, two-to-three-rank restart and both serial/MPI
restart directions produce byte-identical state CSVs. Same-rank restart also
matches budget CSV bytes. Independent inventory/budget checks pass with maximum
normalized residual **3.18097923458357769e-15**. Changed automatic mode and CFL
are rejected with the expected policy diagnostics while protecting output
sentinels; corrupt record, changed cap and unsupported subcycling rejection
checks also pass.

Configure/build, provenance, Python syntax, project-contract and whitespace
checks pass. Compiler warnings remain. No full matrix, selected mechanism,
installed-binary, regridding/subcycling, scaling or physical validation was run.

Logs: `build/mpi-eb-tree-auto-checks/{configure,build,ctest}.log`.
Runtime artifacts: `build/mpi-eb-tree-auto-debug/tests/mpi_eb_tree_auto_3d/`.
