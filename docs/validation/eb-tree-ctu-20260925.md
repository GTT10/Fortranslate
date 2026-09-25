# Extended EB patch geometry and serial CTU tree stepping — 2026-09-25

Implementation `f8e37f9`, parent `19ebfae`, introduces
`build_eb_tree_extended_geometry_3d`. The finest full-domain geometry must be
padded by the requested patch width times the total refinement ratio (default
patch width six). Existing conservative geometry coarsening builds each level;
the builder extracts patch margins, checks their central metrics against the
accepted physical patches, and retains those physical metrics bit for bit.
This handles clipping-coordinate roundoff without replacing physical geometry
or weakening the existing matching tolerance. Outputs publish transactionally.
The builder takes topology and replicated geometry, independent of state ownership.

`advance_reactive_eb_tree_hydro_3d` accepts optional extended patch geometry.
Supplying it selects `pelec` PLM CTU, six-layer peer/ancestor halos,
multidirectional flattening, physical outflow/Mach face positions, viscosity,
centroid fluxes and integrated-flux publication. Existing peer reconciliation,
coarse/fine reflux, restriction and composite boundary/wall accounting follow.
The common interval still uses PCM state donors and unredistributed patch
updates, as the preceding tree hydro path did. Omitted geometry retains PCM
Euler. Covered storage is preserved and all patches commit together.
The algorithm reference remains PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

## Checks

GNU Fortran 13.3.0 / OpenMPI 4.1.6 / Python 3.12.3; existing Ninja Debug
`build/spray-amr-reactive-debug`. No compiler/backend changes.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_hydro_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_(hydro|viscosity)_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final **2/2 passed**: explicit-zero coefficient 14.19 s, coefficient .037
14.26 s; 28.49 s total. Each exercises x/y/z cuts with four- and five-node
three-level trees. Accepted CTU runs start from discontinuous densities at
1000 K, with restriction before stepping. Central cache volume fractions are
exactly equal to physical patches. Deep-node states change, and composite
inventory residuals after boundary/wall accounting remain below the unchanged
3e-14 tolerance; maximum observed CTU residual **5.999e-15**. Incompatible solver
and the cold Riemann fixture reject without changing state or ledgers.
The new CTU fixture begins with zero velocity, so positive-viscosity execution
here does not establish a nonzero artificial-viscosity effect. Existing PCM
and Mach tests remain in the same programs.

### Initial failures and thermodynamic boundary

The first trial reused the old cold pressure-jump fixture and failed: its
Riemann star temperature was **284.4554 K**, below the mechanism's common
minimum **300 K**. A subsequent pressure-doubled fixture passed x/y but
failed z because a reconstructed primitive state exceeded the mechanism range
(observed rho .31, pressure 540000). The final success fixture uses 1000 K
at every active input cell, retaining density discontinuities. The original
cold x/four-node case remains an explicit rejection/rollback test. No EOS
bounds, numerical tolerances or golden data were loosened. Temporary diagnostic
prints were removed. There were no build failures or skips.

The diagnostic runs also raised an overflow flag in the existing shock-sensor
zero-denominator path; the shared sensor's denominator safeguard remains a
separate numerical follow-up. Compiler warnings remain.

Artifacts in `build/eb-tree-ctu-checks/` include `build.log`, `rebuild.log`,
`bounded-build.log`, `LastTest-first.log`, `LastTest-warm.log`, final
`LastTest.log`, and the `debug-*.log` diagnostic builds/runs.
Provenance matched SHA
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`;
Python syntax, project contract and diff checks passed.

## Remaining scope

This is a serial common-step library path. Recursive subcycling, owner-local
MPI stepping, application/restart method selection and coupled tree
reaction/transport dispatch remain open. Other geometry families, nonzero CTU
viscosity effects and wide Mach CTU runs are not newly qualified here. No
application matrix, installed-binary audit, external solver parity or physical
validation was performed.
