# Repeated wall reflections in wide AMR halos — 2026-09-25

Implementation `f266274`, based on `a7e84dc`, removes the single-domain-width restriction
from the shared serial/MPI AMR donor sampler. Ghost coordinates reflect until
they reach an interior donor or an open boundary. Velocity transforms are
composed from the donor outward, preserving the order of distinct moving-wall
velocities. Total energy retains the existing kinetic-energy correction.
The public halo-width bound remains 1 through 128. Default two-layer behavior
and transactional publication are unchanged.

This supplies a prerequisite for wide-filter Dynamic LES at walls; it does
not enable Dynamic wall dispatch or implement an isothermal heat-flux model.

## Focused evidence

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug, existing fixed tests-enabled
configuration. Both changed test targets built successfully.

Final CTest: **2/2 passed**, no failures or skips, total 8.57 seconds:

- `unit_amr_patch_tree_ghost_3d`: 0.09 seconds. A 13-layer halo around a
  four-cell domain matches an independent analytic extension for slip walls,
  no-slip walls with distinct tangential velocities (1 and 3), and a wall
  opposite an extrapolated outlet. All cells compare donor values, momentum
  and total energy exactly. Existing dense-grid and temporal donor oracles
  also pass.
- `unit_mpi_amr_patch_tree_ghost_3d_2`: 8.44 seconds. Seven-layer periodic,
  mixed physical and paired moving-wall halos match serial exactly for PCM,
  MC and minmod reconstruction. Distributed and concentrated ownership,
  empty ranks, temporal donors and collective rollback remain covered.

An earlier test revision also passed 2/2 in 5.05 seconds, before adding the
paired moving-wall MPI fixtures. No tolerance was changed. Provenance,
Python syntax, project contract and whitespace checks passed.
This is boundary-sampling evidence, not LES wall qualification or physical
validation. Selected-mechanism builds and a full matrix were not rerun.

## Commands and artifacts

The final CTest log is copied to `build/amr-wide-wall-checks/LastTest.log`.
MPI uses the previously recorded host-local launcher workaround.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_amr_patch_tree_ghost_3d test_mpi_amr_patch_tree_ghost_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_(mpi_)?amr_patch_tree_ghost_3d(_2)?$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
