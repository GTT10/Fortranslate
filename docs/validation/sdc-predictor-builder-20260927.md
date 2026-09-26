# MPI EB SDC predictor builder — 2026-09-27

Implementation `8143126a`, based on `d360865a`, local `main`, version 0.246.0.
Pinned PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Implementation

`advance_mpi_eb_tree_sdc_3d` accepts optional `predictor_builder` with public
interface `mpi_eb_sdc_predictor_builder`. Each iteration passes the combined
valid-cell predictor rates, interval, iteration number and selected depth to
the caller. The output is a fresh complete three-cell halo, passed to hydro's
existing ownership, shape, finiteness and valid-fluid-cell equality checks.
The builder must include every applicable source family. It may capture
geometry, temporal donors and particle stage data in its own context.

MPI agrees on callback presence and callback success. A failure on any rank
returns before publishing gas, temperature, reaction estimates or inventories.
Callback-owned data remain the caller's transactional responsibility.
Without the callback, hydro retains its normal source fill.

This connects a caller-assembled halo to the level/common-step SDC iteration.
Recursive-driver forwarding and actual physical/virtual/ghost particle source
assembly remain open; this is not an integrated recursive spray application.

## Completed checks

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6, `/usr/bin/gfortran`, `/usr/bin/mpiexec`.
`unit_mpi_eb_tree_sdc_aux_source_3d_2` passed: 1/1, 48.33 s; total 48.38 s.
No build or test failures and no test skips occurred.

The existing three-level reacting/transport/auxiliary-source fixture executes
both distributed and all-on-rank-zero ownership. The builder runs exactly once
per iteration and receives the expected interval/depth. Normally filled halos
retain exact serial state, temperature and reaction-source agreement, with the
existing inventory tolerances unchanged. A builder failure on only the last
rank during iteration two collectively rejects and preserves gas, temperature,
reaction source, all checked source/boundary/wall/chemical inventories and limiter.
Existing three late auxiliary-source rejection cases also pass.

The focused check uses common-step depth -1 and an unmodified normal halo.
Selected-level callback use and modified particle halos were not independently
run here. No full matrix, installed audit, external solver numerical parity or
physical validation was performed. Provenance, Python syntax, project contract
and whitespace checks passed.

## Reproduction

Artifacts: `build/sdc-predictor-builder-checks/` (build, provenance, contract,
CTest summary and copied detailed `LastTest.log`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_full_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_sdc_aux_source_3d_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
