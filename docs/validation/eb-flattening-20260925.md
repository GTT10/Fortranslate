# EB flattening preparation — 2026-09-25

## Scope and source

Implementation `62cabdc`, parent `1bf4461`, adds
`compute_reactive_eb_flattening_3d`. It builds the one-layer donor coefficient
field from three-layer conserved-state, temperature and volume-fraction halos.
Each cell takes the minimum sensor over three directions. A missing seven-cell
pencil or a pencil containing a cut cell retains unity in that direction.
Covered pencils also retain unity to avoid reading undefined thermodynamics;
unused exterior edge/corner output donors retain unity. Failure preserves the
caller output. State halos must include the tangential pencils of face donors,
which can read exterior edge/corner input cells.

Reference: pinned PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`,
`Source/Godunov.H:188-298` (`flatten`, `flatten_eb`) and
`Source/EB.H:189-208` (`is_cut_neighborhood`). The existing regular Fortran
sensor is reused, including its pressure-floor and denominator safeguards;
this is not a new compiled external-oracle parity claim. Geometry in this
API is represented by volume fractions, not arbitrary AMReX cell flags.

## Focused evidence

GNU Fortran 13.3.0 / OpenMPI 4.1.6 / Python 3.12.3; existing Ninja Debug
configuration `build/spray-amr-reactive-debug`. No compiler/backend change.

Commands from the repository root:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_eb_hydro_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^regression_reactive_eb_hydro_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final focused result: **1/1 passed**, 8.30 s test / 8.34 s total.
The regression checks a compressive pressure jump in each axis, each of seven
cut-neighbor positions per axis, smooth/unused donors, poisoned covered
storage and transactional failure for invalid active state. The prepared CTU
update now receives computed flattening for regular and three axis-cut
geometries, with all three redistribution choices; existing inventory and
species closure checks pass unchanged.

The first build failed because the new test passed an extra temperature-guess
argument to primitive-to-conserved conversion. That call was corrected. A
CTest invocation before checking this failed build ran the old executable;
its result is **excluded**. The rebuilt sensor test passed (8.01 s), then the
update integration was added and the final check above passed. No tolerance
or golden-data changes. Compiler warnings remain.

Artifacts: `build/eb-flattening-checks/build.log`, `rebuild.log`,
`integration-build.log`, `LastTest-stale.log` (excluded), `LastTest-first.log`
and `LastTest.log`. Provenance SHA matched
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`;
Python syntax, project contract and diff checks passed.

Application geometry/halo preparation, physical boundary corrections,
input selection and MPI/AMR dispatch remain open. No application matrix,
installed-binary audit, external solver run or physical validation was run.
