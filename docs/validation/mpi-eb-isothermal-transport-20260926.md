# MPI tree isothermal transport — 2026-09-26

Implementation 1725502, parent d44627e, local main, version 0.246.0.
MPI common/recursive transport now accepts optional wall_temperature(6),
validates finite positive values and collectively agrees presence and values
before updating private state. Isothermal faces require conduction and supplied
temperatures. Only domain faces pass thermal labels to patch transport.
MPI halo/flux controls now distinguish isothermal walls, with no-slip reflection
kinematics and thermal flux handled by transport.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Python 3.12.3, Ninja Debug. Existing fixed
tests-enabled build/spray-amr-reactive-debug, no backend changes.
Two direct three-rank MPI runs passed, common and recursive transport, no
failures or skips. The five-patch z-planar fixture has all exterior walls at
2000 K. Distributed and all-on-rank-zero ownership exercise empty ranks.
State, temperature, boundary ledger and limiter match serial exactly.
Mass/species boundary flux is zero and heat input is positive. Maximum MPI
inventory residuals: common 2.3158526836485040e-15, recursive
2.1262118844211792e-15, both below the existing 3e-14 bound.
A rank-local 2001 K override rejects and preserves public state/ledgers; existing
invalid dt/state/transport-control rollback cases remain active.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_transport_3d test_mpi_eb_tree_subcycle_transport_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_transport_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none /usr/bin/mpiexec -n 3 build/spray-amr-reactive-debug/tests/test_mpi_eb_tree_transport_3d isothermal
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none /usr/bin/mpiexec -n 3 build/spray-amr-reactive-debug/tests/test_mpi_eb_tree_subcycle_transport_3d isothermal
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The second build includes corrected placement of the fixture's temperature
reset before the second ownership layout. No failed runtime attempts occurred.
Matching CTest entries were added; the reported runs used mpiexec directly.
Provenance, syntax, project contract and whitespace passed. Compiler warnings
remain; no tolerance or golden-data changes.
Artifacts: build/mpi-eb-isothermal-transport-checks/{build.log,build-final.log,
common.log,recursive.log}.
Mechanism SHA: 0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4.
Pinned PeleC: bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8.
Pinned PelePhysics: 6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b.

## Limits

Thermal timestep bounds, full reacting temperature forwarding and application
input/restart remain open. These runs use a fixed small timestep; existing
nonthermal timestep checks do not qualify an isothermal bound. No full matrix,
installed audit, external solver parity or physical validation ran.
Complete PeleC parity remains unproven.
