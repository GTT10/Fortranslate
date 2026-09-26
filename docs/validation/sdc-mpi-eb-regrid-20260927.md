# MPI EB tree SDC regridding — 2026-09-27

Implementation `6de6c040`, based on `9d2a86ef` on local `main`.
PeleC pin: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

`remap_mpi_eb_tree_3d` now accepts an optional owner-local reaction-source field.
Presence is agreed across ranks. A private source copy is restricted from the
deepest descendants using fluid-volume weights, then follows the gas donor map
through same-level overlap or ancestor PCM interpolation. Solid target values
are zero. Final source restriction precedes publication of source, gas,
temperature and geometry. The signed field never passes through EOS recovery,
clipping or positivity enforcement. All stored input source values must be finite.

The reusable `average_down_mpi_eb_tree_field_3d` computes restriction on child
owners and sends coarse signed values to parent owners; there is no root gather.
MPI EB application regridding now passes the SDC source field and refreshes its
checkpoint copy after success. The prior SDC regrid/adaptation input gate is
removed; recursive SDC stepping remains unavailable.

## Completed focused checks

GNU Fortran 13.3.0, `/usr/bin/mpifort` and `/usr/bin/mpiexec`
(OpenMPI/OpenRTE 4.1.6), Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`. Existing compiler warnings remain.

- Mechanism provenance, relevant builds, Python syntax, project contract and
  `git diff --check` passed.
- CTest `unit_mpi_eb_tree_remap_3d_2`: 1/1 passed, 1.34 s test / 1.39 s total.
- The existing MPI remap fixture now carries negative, spatially varying source
  estimates through patch movement, complete coarsening and three-level
  refinement while changing ownership. Gas, temperature and sources agree with
  serial exactly; sources remain owner-local. Source inventory is conserved
  within the asserted normalized bound 3e-14.
- Invalid gas and changed geometry reject without changing gas, temperature,
  geometry or reaction estimates.
- Four application runs passed: serial and two-rank continuous, two-rank stop
  after the first regrid, and one-rank restart. The SDC fixture runs three steps
  to 3e-10 s, coarsens 5 patches to 1 at 1e-10 s and refines to 5 at 2e-10 s.
- Serial/MPI state CSVs are byte-identical. Two-rank continuous and one-rank
  resumed state and budget CSVs are byte-identical. The MPI composite residual
  is `1.7210356769130377e-16`; the maximum normalized residual across all checked
  serial/MPI budget CSVs is `8.389084975536366e-15` (asserted bound 5e-10).

Commands from repository root:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_remap_3d pelef_mpi_reactive_eb_tree_3d pelef_reactive_eb_tree_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_remap_3d_2$' --output-on-failure
python3 build/sdc-mpi-eb-regrid-checks/run.py
python3 -m py_compile tools/check_project_contract.py build/sdc-mpi-eb-regrid-checks/run.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/sdc-mpi-eb-regrid-checks/`, including the local runner, generated
inputs, logs, checkpoint, CSVs and summary. The runner uses the configured
OpenMPI launcher, clears DISPLAY/LD_LIBRARY_PATH/LIBRARY_PATH, and sets the same
HWLOC and oversubscription environment as the unit command.
No automatic-tagging trajectory, simultaneous measured rebalance/regrid,
selected build, full configuration matrix, installed-binary audit, recursive
SDC integration or physical validation was run.
