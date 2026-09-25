# EB tree adiabatic wall transport — 2026-09-26

Implementation e17d0ea, parent 23104f5, local main, version 0.246.0.
The patch Euler transport API now forwards optional physical-face labels to
existing wall flux construction. Serial/MPI common and recursive tree paths
mark only domain faces as slip/no-slip walls, using the integer domain layout;
internal patch interfaces remain ordinary donor interfaces. Periodic labels
remain handled by tree halos and synchronization, not local patch wrapping.

Existing wall flux handling enforces zero mass/species and adiabatic energy
flux, while retaining appropriate viscous momentum stress. This avoids treating
reflected ghost states as unconstrained open diffusive faces.

## Completed focused checks

GNU Fortran 13.3.0, Python 3.12.3, Ninja Debug and existing OpenMPI wrappers.
Fixed tests-enabled build/spray-amr-reactive-debug; no backend changes.
Serial common/recursive transport tests passed **2/2**, 1.85 and 16.33 seconds
(CTest total 18.22 seconds), no failures or skips. Three plane orientations,
four/five-patch trees and all exterior no-slip walls are exercised. Boundary
mass, species and energy fluxes are exactly zero. Component inventory residuals
remain below 3e-14; maximum printed wall residual is 5.2293019464126058e-15.
Existing periodic, outflow and failure checks remain in these tests.
MPI common/recursive source forwarding was compiled through the MPI gas
application; MPI wall transport runtime was not checked in this increment.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_transport_3d test_eb_tree_subcycle_transport_3d pelef_mpi_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_transport_3d test_eb_tree_subcycle_transport_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_(transport|subcycle_transport)_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The second build incorporates the final closed-wall flux assertions.
Provenance, syntax, project contract and whitespace passed. Existing compiler
warnings remain. No tolerance or golden-data changes.
Artifacts: build/eb-tree-wall-transport-checks/build.log, build-tests.log and
LastTest.log. Mechanism SHA:
0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4.
Pinned PeleC: bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8.
Pinned PelePhysics: 6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b.

## Limits

MPI runtime, slip-wall transport stress qualification, full reacting wall split,
application wall controls, general reflected CTU geometry and isothermal tree
walls remain separate work. No full matrix, installed audit, external solver
parity or physical wall-flow validation ran. Complete PeleC parity is unproven.
