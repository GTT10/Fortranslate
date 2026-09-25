# Reacting MPI EB tree wall splits — 2026-09-26

Fixture 5619881, parent 78c6ebb, local main, version 0.246.0.
This adds wall modes to the existing public full-split test and reruns serial
and MPI reaction/transport/hydro composition. No solver change was needed.
It completes focused runtime evidence for the MPI wall transport forwarding
introduced in e17d0ea, together with the existing wall hydro/halo connections.

## Completed evidence

GNU Fortran 13.3.0, Open MPI 4.1.6, Python 3.12.3, Ninja Debug, existing
fixed tests-enabled build/spray-amr-reactive-debug. No backend changes.
Three-rank recursive/common CTU tests passed **2/2**, 111.88 and 89.57 seconds
(CTest total 201.50 seconds), no failures or skips.

The five-patch, three-level z-planar fixture uses x-slip, y-no-slip and z-outflow
faces. Explicit chemistry uses distributed owners; implicit chemistry uses the
last rank only, leaving two ranks empty. Both include molecular transport and
hydro. MPI state, temperature, boundary ledger and limiter equal the serial
results exactly. Embedded-wall/chemical ledgers and inventory satisfy existing
3e-14 bounds. Printed MPI inventory residuals are about 4.30e-15. Serial/MPI
boundary-aware timestep estimates satisfy the existing 1e-14 relative check.
A rank-local wall-label mismatch after private reaction work rejects without
changing public state or ledgers. Existing other failure cases remain active.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_subcycle_full_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_ctu-wall(-common)?_full_3d_3$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, project contract and whitespace passed. Compiler warnings
remain. No tolerances or golden data changed.
Artifacts: build/mpi-eb-wall-full-checks/build.log and LastTest.log.
Mechanism SHA: 0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4.
Pinned PeleC: bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8.
Pinned PelePhysics: 6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b.

## Limits

Walls are tangential to the planar EB. General reflected CTU geometry,
application wall input/restart, spray wall routing and isothermal tree walls
remain separate work. No full matrix, installed audit, external solver parity
or physical wall-flow validation ran. Complete PeleC parity remains unproven.
