# MPI EB tree Mach inlet stepping — 2026-09-25

Implementation `44bfe3e`, base `cd10be4`.
The owner-local EB halo API now accepts the same optional six-face
Mach, temperature and pressure controls as the serial halo. It validates finite
values and active-inlet ranges collectively, compares all 18 effective control
values across ranks, then constructs physical inlet ghosts from the fetched
interior donor. Covered faces, peer/coarse donors and temporal interpolation
retain their previous treatment. Failed state construction does not publish halos.

Common and recursive MPI hydro/transport plus the complete reactive split forward
the controls. Active inlet faces permit molecular flux, while inactive physical
faces retain closed transport. Existing local candidate-state commit, reflux,
boundary/wall/chemical ledgers and empty-owner participation remain.
No gathered global state or new persistent checkpoint field is introduced.

The shared Mach state kernel retains the pinned nozzle equations described in
[the state evidence](mach-inlet-state-20260925.md). This increment tests MPI/serial
agreement, not a new comparison with the reference C++ solver.

## Focused verification

GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpiexec`), existing Ninja Debug
fixed-mechanism build. The five-patch plane-cut tree has two refinement depths.
Low-x Mach 0.2 extrapolates pressure; high-y Mach 0.1 prescribes 120000 Pa;
both use 400 K. Two ranks advance 1e-8 s with chemistry and all molecular
transport terms enabled, in both common and recursive modes. Explicit chemistry
uses round-robin patch ownership; implicit chemistry places all patches on the
last rank, exercising a rank with no owned cells.

Tests require exact serial/MPI conserved states, temperatures, physical-boundary
ledger and limiter. Chemical/wall ledger parity and total conservation retain
3e-14 tolerance. A one-rank pressure discrepancy and one-rank NaN temperature
must reject collectively with state and diagnostic rollback. Existing late
hydro failure and inconsistent chemistry/subcycling flags are also exercised.
No tolerances changed. The Mach variants run the fully coupled mode only;
the pre-existing inactive-inlet test retains its original mode combinations.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_subcycle_full_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_mach(-common)?_full_3d_2$' --output-on-failure --output-log build/eb-mach-mpi-checks/ctest.log
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**: recursive **42.87 s**, common **34.45 s**, total **77.37 s**.
Build, provenance, Python syntax, project contract and whitespace passed.
No failures or skips; compiler warnings remain. Artifacts:
`build/eb-mach-mpi-checks/{build.log,ctest.log,LastTest.log}`.

Inlet-aware timestep bounds and application input/checkpoint wiring are still
pending. No selected-mechanism execution, full matrix, large-tree scaling,
changed-rank restart, external solver run or physical nozzle validation was
performed for this library increment.
