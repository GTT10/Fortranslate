# MPI EB SDC diffusion callback — 2026-09-27

Implementation `36e4678f`, based on `8cf2332b`, local main, version 0.246.0.

## Change

Common and recursive MPI EB SDC accept `evaluate_diffusion`. The callback
receives owner-local gas state/temperature, the local interval and selected
level, plus parent old/new state/temperature donors and stage fractions.
It returns rates per fluid volume/time, interval-integrated fluxes per open
face area, and a limiter value. Normal face indices must start at zero, matching
`eb_integrated_flux_3d`; aperture factors remain the synchronization operator's
responsibility.

Callback rates/fluxes add to molecular transport before existing redistribution.
Combined old/new fluxes participate in final synchronization even when native
transport is disabled. Presence, status, owner-local layout, active rate
finiteness, flux bounds/shapes/finiteness and limiter range are checked
collectively. The limiter is reduced by the global minimum. Recursive calls
forward the callback through all auxiliary-source/predictor combinations.
Callback physical state, if any, must remain private until outer-step success.

## Independently completed checks

Artifacts: `build/eb-sdc-diffusion-checks/` (`provenance.log`, `build.log`,
`build-final.log`, `ctest.log`, `ctest-final.log`, `LastTest-initial.log`,
`LastTest-final.log`). Existing fixed-mechanism Ninja Debug build, GNU Fortran
13.3.0, Open MPI 4.1.6 (`/usr/bin/mpiexec`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_full_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_sdc_diffusion_3d_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final **1/1 passed**, 91.95 s (CTest total 92.00 s), no failures/skips.
Build, mechanism provenance, Python syntax, contract and whitespace checks
passed. Existing real-comparison warnings and executable-stack linker warning
remain. Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

The two-rank, five-node EB tree test routes the existing molecular transport
operator through the callback with native transport disabled. Common-step gas,
temperature and reaction estimates match the serial native-transport reference
exactly; inventories meet the existing 2e-12 comparison. Composite normalized
residual is 4.2614898476812615e-15 in both distributed and rank-0-only layouts.
Recursive callback state, temperature, reaction estimates and all diagnostics
match native recursive MPI transport exactly. Callback failure on the last rank
at its second evaluation rejects the whole recursive step and preserves those
outputs, including when that rank owns no cells.

An initial compile rejected a repeated host-associated procedure `present`
query; callback presence is now retained in a host logical. The first test
failed because its adapter passed `target_depth=-1` explicitly to a transport
API that requires omission for common steps. The adapter now omits that option
for common steps. Initial failure logs remain separate from final evidence.

## Limits

This checks the gas integration interface using molecular transport, not LES
or particle coupling. An actual LES adapter and particle-wrapper/application
forwarding remain open. No cut-cell LES model was added. Combined nonzero
native/callback diffusion, callback redistribution cases, malformed outputs,
selected mechanisms, installed binaries, configuration matrices, external
PeleC parity and physical validation were not independently tested here.
