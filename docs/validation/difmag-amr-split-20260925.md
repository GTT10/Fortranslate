# AMR coupled artificial-viscosity propagation — 2026-09-25

## Change

Implementation `80c7a85`, parent `3092723`. Serial/MPI patch-tree
reaction/transport/hydro splitting forwards optional `difmag` through common or
recursive hydro. The reaction/hydro Strang wrappers also forward it. The MPI
AMR spray coupled wrapper passes it into the gas split. Transport and chemistry
do not consume the artificial-viscosity coefficient themselves.

Validation and collective coefficient agreement remain in hydro. Earlier stage
updates are private candidates, so rejection preserves the public gas state,
temperature and diagnostics. This connection does not yet enable coefficient
input in AMR applications or add restart identity.

## Verification

GNU Fortran 13.3.0, OpenMPI 4.1.6, Ninja Debug, Python 3.12.3; fixed-mechanism
tests-enabled `build/spray-amr-reactive-debug`.

**2/2 passed**, common fixture 73.24 s and recursive fixture 128.19 s,
201.48 s total elapsed. Build, provenance, Python syntax, project contract
and whitespace checks passed. No failures or skips occurred in this increment.
Compiler warnings remain.

The new PPM variants of existing split fixtures use `pelec` and coefficient
0.037. Conserved state, temperature, boundary ledger, reflux maximum and
transport limiter remain exactly equal to serial. Fixtures include periodic,
inflow/outflow, moving isothermal walls, pressure outlets and owner maps with
an empty rank. The common fixture additionally includes Mach inlets and a
recursive Mach case. Explicit and implicit chemistry and all molecular
transport terms run in enabled cases; a disabled-operator case remains.
Rank-dependent coefficient 0.037/0.05 rejects the full composition, preserving
state and diagnostics after preceding chemistry and transport. Existing bad
hydro-solver and inconsistent-operator guards also run without relaxed checks.

The AMR spray executable builds with the forwarded optional argument. This is
compile evidence for spray, not a new particle runtime check. The standalone
reaction/hydro Strang wrappers were compiled, not independently exercised with
a custom coefficient. No selected-mechanism build, application/restart test,
full matrix, installed audit, external solver comparison or physical validation
was performed in this increment.

## Commands and artifacts

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_split_3d test_mpi_amr_patch_tree_subcycle_split_3d pelef_mpi_amr_spray_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_amr_patch_tree_ppm_(split|subcycle_split)_3d_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Logs: `build/difmag-amr-split-checks/build.log`, `test.log`, `LastTest.log`.
