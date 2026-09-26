# MPI EB spray common-step SDC application — 2026-09-27

Implementation `5cd1a8c4`. Final CTest **1/1 passed**, 135.40 s total
(test 135.35 s). Provenance, build, Python syntax, project contract and
`git diff --check` passed. No solver/build failures occurred.
Base: `f6ce6169` on local `main`; version remains 0.246.0.
PeleC pin `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Change

The MPI EB spray application now dispatches `time_integrator='sdc'` to the
common-step particle/gas wrapper. It passes reaction history, gas physical
boundary data, redistribution controls, particle IDs and split counts.
Gas acceleration is `external_forcing + acceleration`; particle acceleration
is the spray `acceleration`. CTU/PeleC hydro and explicit faces remain required.
Film coupling and recursive particle SDC still reject at input validation.

Coupled restart saves total external work in its prefix even when the spray
acceleration is zero. The embedded gas snapshot uses a zero-exchange reference
and zero forcing history but retains its reaction source. Restore publishes
that source alongside gas/particles and restores the coupled external inventory.
Existing gas checkpoint identity and spray context bind integrator controls and
acceleration. Default split behavior is retained.

## Evidence boundary

GNU Fortran 13.3.0, Ninja Debug, fixed-mechanism tests-enabled build
`build/spray-amr-reactive-debug`; OpenMPI 4.1.6 through `/usr/bin/mpifort`
and `/usr/bin/mpiexec`. The existing executable-stack linker warning remains.

The checker uses three 1e-8 s steps, two AMR levels, chemistry and transport,
water parcels, EB reflection, gas-only forcing `[100,-200,300]` and common
acceleration `[1000,-2000,3000]` m/s². It compares two/four-rank gas and sorted
particle outputs; stops after step two and continues with four ranks; compares
restarted gas and particle CSV values exactly and all budget components,
including `force`, within 3e-12 relative. Missing manifests, attempted archive
overwrite, mismatched spray context/hydro and truncated particle payload are
also rejected. These are focused application/restart checks, not physical
spray validation or complete PeleC parity. Selected-mechanism execution,
SDC spray regridding/rebalancing, films and recursive stepping were not tested.

The first checker invocation completed continuous and restarted calculations,
but its new budget assertion used `forcing` instead of the spray CSV column
`force`, raising a Python KeyError. The checker was corrected without changing
the solver or tolerances. Independent inspection of those completed artifacts
found the full restarted budget CSV exactly equal to the continuous four-rank
budget; each run reported normalized composite error 1.2834732843575932e-17.
The final CTest rerun includes the remaining rejection checks.

Artifacts: `build/spray-sdc-app-checks/` (first invocation and build logs) and
`build/spray-amr-reactive-debug/spray-sdc-restart/` (final CTest cases).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc_restart$' --output-on-failure
python3 -m py_compile tests/application/check_mpi_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
