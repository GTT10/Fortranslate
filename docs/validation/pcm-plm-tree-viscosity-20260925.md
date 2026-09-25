# PCM/PLM tree hydro artificial viscosity — 2026-09-25

## Source and scope

Code `24a4fbb`, parent `3ce8b6e`. Serial and MPI patch-tree common/level hydro
now accept positive explicit PCM/PLM coefficients and apply the existing
regular-cell correction after reconstruction fluxes, before cell updates,
time-integrated flux storage and reflux. Existing recursive drivers forward
the coefficient. Physical no-slip masks come from the existing domain-touch
helper; internal patch faces are not treated as walls. MPI retains collective
finite/nonnegative coefficient validation and exact control agreement.

Omitted/zero PCM/PLM coefficients retain their old policy. Application input,
default/restart policy, older two-level AMR, regular MPI spray and EB remain
outstanding. Coupled chemistry/transport wrappers already forward coefficients,
but this increment tests hydro only. This is not full PeleC solver parity or
physical validation.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Ninja Debug, Python 3.12.3; reused
`build/spray-amr-reactive-debug`. Focused CTest **2/2 passed**: common hydro
6.40 s; recursive subcycling 20.02 s; total elapsed 26.46 s. No failures or
skips. Compiler real-comparison warnings remain.

Both two-rank fixtures use the `pelec` solver and coefficient 0.037 across
eight existing boundary/ownership modes alternating PCM and PLM. Serial and
MPI conserved states, temperatures, boundary budgets and reflux diagnostics
match exactly. Recursive level counts match 1/2/6. A coefficient mismatch
(0.037 versus 0.05) and a negative coefficient on one rank reject collectively
without changing state; existing failure/rollback checks also pass. The
fixtures cover periodic, slip-wall, outflow, pressure-outlet and inflow cases;
they do not separately exercise tree no-slip viscosity in this increment.

Mechanism provenance, Python syntax, project contract and whitespace checks
passed. No tolerances were loosened.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_hydro_3d test_mpi_amr_patch_tree_subcycle_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_amr_tree_viscosity_(hydro|subcycle)_3d_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/pcm-plm-tree-viscosity-checks/build.log` and copied
`LastTest.log`.
