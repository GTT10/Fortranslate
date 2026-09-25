# MPI periodic EB tree hydro — 2026-09-26

Source `c5a2981`, parent `e671ddc`, local main. Common and recursive MPI hydro
forward `boundary_faces` into flux synchronization and no longer reject
periodic labels up front. Recursive stepping also forwards them into the
peers-only correction before publishing provisional ancestor endpoints.
Existing halo and synchronization validation enforce collective agreement.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Python 3.12.3, Ninja Debug; existing
tests-enabled build, no toolchain changes.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_hydro_3d test_mpi_eb_tree_subcycle_hydro_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_periodic_(hydro|subcycle_hydro)_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, common 6.89 s, recursive 18.82 s, total 25.76 s.
No failures or skips. Both use three ranks, five patches, three levels,
PLM/CTU with nonzero artificial viscosity, x/y periodic and z outflow faces.
Distributed owners and all patches on the last rank (two idle ranks) match
serial state, temperature and boundary inventory arrays exactly. Wall ledger
and composite conservation satisfy the unchanged 3e-14 tolerance.
Existing malformed solver, timestep, state, geometry and method checks retain
rollback assertions. A new rank-zero boundary-label disagreement also rejects
without changing owner-local state, temperature or diagnostic outputs.

Provenance, Python syntax, project contract and whitespace checks passed.
Artifacts: `build/mpi-eb-periodic-hydro-checks/build.log` and `LastTest.log`.
Mechanism SHA:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

MPI periodic molecular transport, reacting splits and application input remain
pending. This check does not establish full configuration coverage, installed
binary qualification, external solver parity or physical validation.
