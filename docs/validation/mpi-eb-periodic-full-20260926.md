# MPI periodic EB tree transport and reacting splits — 2026-09-26

Source `54f988d`, parent `21d8b9f`, local main. Common and recursive molecular
transport now forward periodic labels to final flux synchronization. Recursive
transport also passes them into provisional endpoint correction before child
steps. Existing full reacting split forwarding therefore connects periodic
reaction/transport/hydro composition without changing its transaction model.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Python 3.12.3, Ninja Debug; existing
tests-enabled build, no compiler/backend changes.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_subcycle_full_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_ctu-periodic(-common)?_full_3d_3$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, recursive 111.62 s, common 87.88 s, total 199.55 s.
No failures or skips. Three ranks, five patches and three levels, PLM/CTU,
x/y periodic and z outflow faces. Both tests cover explicit/implicit chemistry
and all four chemistry/transport enablement combinations. Explicit cases
distribute patches; implicit cases put all patches on the last rank with two
idle ranks. Owner-local state, temperature, boundary inventories and limiter
match serial exactly. Wall/chemical inventories and total conservation satisfy
the unchanged 3e-14 bounds. Donor timestep bounds agree to the existing 1e-14
relative tolerance. Unpaired labels on rank zero reject after private reaction
without publishing state or diagnostic changes; existing split-failure checks
remain active.

Provenance, Python syntax, project contract and whitespace checks passed.
Artifacts: `build/mpi-eb-periodic-full-checks/build.log` and `LastTest.log`.
Mechanism SHA:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

Tree application periodic input, restart/regrid qualification and MPI spray
periodic coupling remain pending. No full matrix, installed audit, external
solver parity or physical validation was run.
