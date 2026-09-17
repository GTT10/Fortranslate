# Owner-local MPI EB tree hydro — 2026-09-18

Implementation: `a064423760012f8ffe364b8c13c08dfc6e80fe25`, based on
`7feb668`. New modules `mpi_eb_tree_flux_3d_mod` and
`mpi_eb_tree_hydro_3d_mod` are part of `pelef_mpi_support`.

## Implemented behavior

`synchronize_mpi_eb_tree_flux_3d` retains serial lexicographic peer ownership
and finest-to-coarsest traversal. It first reconciles raw same-level transfers,
then corrects exposed cell states from effective-minus-raw flux divergence.
Fine fluxes are restricted by open area, copied into their parent and mirrored
into adjacent parent-level patches, including one-sided deeper refinement.
Optional minimum/maximum levels are agreed across ranks. Final volume-weighted
averaging precedes publication of U, T and flux together.

States, temperatures and persistent fluxes exist only on patch owners. A
source-owner/receiver exchange temporarily transfers one complete patch flux
at a time; other ranks allocate no such payload. This avoids a root gather,
but full-patch flux messages and sequential pair traversal are not a scaling
optimization. Geometry is replicated and collectively validated.

`advance_mpi_eb_tree_hydro_3d` connects private average-down, sparse spatial
face exchange, owner-local Euler updates, flux synchronization, and separate
physical-boundary/slip-wall inventories. Boundary flux comes from the effective
root transfers; wall impulses include only exposed cells. Timestep and solver
must agree across ranks. Any collective rejection preserves both field
containers and both ledgers. Physical boundaries remain outflow. Temporal
subcycling, molecular transport, reaction composition, selected-target wiring
and an MPI EB tree application remain separate work.

## Completed focused verification

Fresh GNU Fortran 13.3.0 Debug/Ninja build; compiler `/usr/bin/gfortran`, MPI
wrapper `/usr/bin/mpifort`, launcher `/usr/bin/mpiexec`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/mpi-eb-tree-flux-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=ON
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/mpi-eb-tree-flux-debug --target test_mpi_eb_tree_flux_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/mpi-eb-tree-flux-debug --target test_mpi_eb_tree_flux_3d test_mpi_eb_tree_hydro_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/mpi-eb-tree-flux-debug -R '^unit_mpi_eb_tree_(flux|hydro)_3d_[23]$' --output-on-failure --output-log build/mpi-eb-tree-flux-checks/ctest.log
python3 -m py_compile tools/check_project_contract.py tools/check_mechanism_provenance.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**4/4 tests passed**, no failures/skips, total **5.02 s**. Flux tests took
0.78/0.81 s; hydro tests took 1.71/1.69 s for two/three ranks respectively.

Both suites cover x/y/z cut planes with four/five-node, three-level trees,
including peers with different parents and one-sided deep refinement. They
exercise distributed owners and all patches on one rank, leaving others empty.

- Manufactured flux synchronization matches serial U/T and every effective
  flux value exactly. The level-1-through-2 scope also matches serial. NaN
  flux, finite nonphysical correction and rank-inconsistent scope fail
  collectively with unchanged U/T and bitwise unchanged flux payloads.
- Common hydro U/T and the boundary ledger match serial exactly. Wall sums
  agree within 3e-14. The **24** independently integrated MPI hydro inventories
  have maximum normalized boundary/wall-adjusted residual
  **5.6005960242831442e-15**, below the unchanged 3e-14 bound. Invalid solver,
  rank-inconsistent timestep and invalid deep state preserve fields/ledgers.

Configure/build succeeded; compiler warnings remain and this is not claimed
as warning-free. Provenance, Python syntax, project contract and whitespace
checks pass. No full configuration matrix, installed audit, selected mechanism,
coupled application, scaling or physical validation was run.

Logs: `build/mpi-eb-tree-flux-checks/{configure,build,ctest}.log`.
Detailed output: `build/mpi-eb-tree-flux-debug/Testing/Temporary/LastTest.log`.
