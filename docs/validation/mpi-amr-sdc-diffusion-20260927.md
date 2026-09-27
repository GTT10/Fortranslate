# MPI AMR SDC diffusion callbacks — 2026-09-27

Implementation `2c63bb8b`, based on `0731ed68`, local main, version 0.246.0.

## Change

The MPI common and recursive SDC drivers accept the serial diffusion callback
interface. Callbacks receive owner-local state/temperature, the selected level,
parent temporal endpoints and stage fractions. They return instantaneous rates
and face flux densities. The driver integrates fluxes over each local interval,
adds them to molecular transport when enabled, and uses the existing old/new
averaging, reflux and recursive synchronization paths.

Ranks must agree on callback presence. Callback failures and invalid owner-local
layouts, flux shapes or nonfinite active rates/fluxes reject the step collectively.
Updates remain private until success. A callback that performs MPI collectives
must itself call those collectives consistently on all ranks, including idle ranks.

## Independently completed checks

Existing Ninja Debug fixed-mechanism build, GNU Fortran 13.3.0,
`/usr/bin/mpiexec` Open MPI 4.1.6. Artifacts:
`build/mpi-amr-sdc-diffusion-checks/` (`provenance.log`, `build-final.log`,
`ctest-final.log`, `LastTest-final.log`; initial build/test logs also retained).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_split_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_amr_patch_tree_sdc_3d_4$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final focused test: **1/1 passed**, 25.83 s (CTest total 25.88 s), no failures
or skips. Initial run passed in 26.05 s; the final run followed addition of
explicit finite-rate validation. Build succeeded with existing warning classes
(real comparisons, unused callback arguments and executable-stack linker warning).
Provenance, Python syntax, project contract and whitespace checks passed.
Pinned H2O2 source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

The four-rank test retains existing coupled chemistry/transport checks and adds
an analytic momentum-flux callback over three levels (ratios 2 and 3), with two
idle ranks. Common steps with molecular transport and recursive steps without
it match serial state, temperature and inventories exactly; recursive counts
also match. Failure on rank 0 at the second callback evaluation rejects a
recursive step collectively and preserves state, temperature and inventories.

## Limits

The analytic callback checks integration plumbing, not an LES discretization.
MPI LES instantaneous RHS/wrapper and application routing remain open. Callback
presence mismatch and malformed/nonfinite callback outputs have code guards but
were not separately injected in this focused test. No configuration matrix,
selected-mechanism audit, installed-binary audit, external PeleC parity or
physical validation was run.
