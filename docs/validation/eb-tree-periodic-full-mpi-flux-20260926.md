# Periodic reacting tree splits and MPI flux synchronization — 2026-09-26

Serial test source `74d3741`, parent `10c5f43`; MPI implementation/test source
`593df62`, parent `74d3741`, local main. Serial split tests add x/y periodic
faces with z outflow. MPI synchronization adds collectively agreed boundary
labels and matching root metric validation. Periodic peer transfers and
restricted fine transfers follow serial ordering, using existing owner-to-owner
payload exchange. Provisional peers-only synchronization also accepts periodic
faces; no MPI time-update caller enables periodic dispatch yet.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Python 3.12.3, existing Ninja Debug
tests-enabled build; versions checked live.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_subcycle_full_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_ctu_periodic_(full|common)_3d$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_flux_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none OMPI_MCA_rmaps_base_oversubscribe=1 /usr/bin/mpiexec -n 3 build/spray-amr-reactive-debug/tests/test_mpi_eb_tree_flux_3d periodic
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Serial CTest: **2/2 passed**, recursive 151.28 s, common 102.26 s,
253.58 s total; no failures or skips. Both exercise explicit/implicit chemistry
and all four chemistry/transport enablement combinations on five patches and
three levels. Composite inventories meet the existing 3e-14 bound; the
independent public-operator split matches all state, temperature and diagnostic
arrays exactly. Invalid hydro and unpaired periodic labels preserve outputs
after private reaction/transport work.

MPI executable: **3-rank PASS**, exit 0. Three plane normals, four/five-patch
trees, distributed owners, all patches on rank zero, and level-scoped
synchronization match serial states, temperatures and effective flux arrays
exactly. Nonfinite flux, nonphysical correction, rank-disagreeing level ranges
and rank-disagreeing periodic labels preserve outputs. Two-/three-rank CTest
entries were added; only the explicit three-rank command above was run here.
The MPI test initially failed to build due to two long source lines; wrapping
them fixed the compile error. No numerical tolerances were changed.

Provenance, Python syntax, project contract and whitespace checks passed.
Artifacts: `build/eb-tree-periodic-full-checks/{build.log,LastTest.log}` and
`build/mpi-eb-periodic-flux-checks/{build.log,initial-build-failure.log,mpi3.log}`.
Mechanism SHA:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

MPI periodic hydro/transport/reacting time updates, peers-only runtime
qualification and application input remain pending. No full matrix,
installed audit, external solver parity or physical validation was run.
