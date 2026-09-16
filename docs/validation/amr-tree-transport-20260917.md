# Regular 3D tree molecular transport — 2026-09-17

Source parent: `2b5f46a293b3107f1f2f57f7ee13ce7ad09b651e` on `main`, plus this change.
PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Transport module SHA256:
`28f38d9984c7a27130192e142aa1c7390afb67c0794d739467d6aacf19cbb6c8`.
Test SHA256:
`451f1fac0caebd4bd09f27e439d74ee7fcdf2fd65901387be79b19ed0a1dddf9`.

## Scope

Common-timestep SSPRK2 transport uses existing viscous, thermal, species and
barodiffusion fluxes, tree viscous halos and nested flux replacement. Global
node geometry selects physical wall faces. Reservoir, tangential moving-wall
and isothermal-wall inputs are forwarded. State/temperature and all returned
diagnostics commit transactionally. The timestep helper includes every node's
local diffusive bound and, when domain/wall arguments are supplied, touched
isothermal-wall bounds.

The existing species limiter is reused; shared-face disagreement or nonphysical
refluxed state rejects the step. No automatic retry, new coarse/fine positivity
limiter, subcycling, PLM coarse/fine donors, combined reaction/transport/hydro
driver, tree MPI/application/restart or physical qualification is claimed.

## Completed checks

GNU Fortran 13.3 Debug with bounds and invalid/zero/overflow checks, OpenMPI 4.1.6;
reused `build/spray-mpi-migration-gnu-debug`, MPI/spray enabled.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target test_amr_patch_tree_transport_3d -j 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^unit_amr_patch_tree_transport_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Final target build and **1/1 focused CTest passed, 0 skips/failures, 10.02 s**.
The test covers root-only comparison against the existing regular-grid transport
update and a four-level/six-node branched tree with adjacent cross-parent boxes.
Boundary modes are periodic, outflow, moving adiabatic no-slip and moving
isothermal walls. Checks include all-node diffusion and wall timestep oracles,
conserved boundary-ledger balance, species closure, positive temperature and
exact state/temperature/diagnostic rollback after a rejected large timestep.

Initial builds caught line-length errors in the new source/test, corrected
before execution. The initial transport test passed; after adding wall timestep
support, the final test above passed again. No tolerances or golden hashes changed.
Provenance and project contract passed; Python AST parsing of the two checker
scripts passed (2/2), as did whitespace checks. Logs:
`build/amr-tree-transport-checks/build.log` and `ctest.log`.
