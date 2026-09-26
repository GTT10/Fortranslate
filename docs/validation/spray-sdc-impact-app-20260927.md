# MPI EB spray SDC film application and restart — 2026-09-27

Implementation `2609cfb5`, based on `d6b2ebc9`, local `main`, version 0.246.0.
PeleC reference `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Change

The MPI EB spray application now forwards film state, wall controls and the
impact selector to the common-step SDC wrapper when `wall_impacts=.true.`.
Empirical-model rejection remains an explicit failure. Existing particle/film
checkpoint V2 carries the resulting films and splash child IDs; the SDC gas
snapshot retains reaction history. Recursive particle SDC remains rejected.

The first application run exposed missing reaction-source remapping in the
spray-specific regrid transaction. After coarsening to a root-only tree, SDC
rejected stale reaction-history fields and exhausted timestep retries. The
spray remapper now takes optional reaction history, remaps a private copy through
the existing gas operator, and publishes it only after particle and film
migration succeed. Application regridding also updates the checkpoint record.
No tolerance, model bound or retry limit was loosened.

## Focused scope

GNU Fortran 13.3.0, Ninja Debug, fixed H2O2 tests-enabled build
`build/spray-amr-reactive-debug`; OpenMPI 4.1.6 through `/usr/bin/mpifort`
and `/usr/bin/mpiexec`. Existing executable-stack linker warning remains.
The checker runs deposition and splash cases, each with four 5e-9 s steps,
chemistry, molecular transport, two SDC iterations, three scheduled regrids
(including root-only coarsening) and requested work rebalancing. It compares
four-rank continuous execution against a two-rank checkpoint resumed on four
ranks, covering gas, sorted particles, films and every budget component.
Two films must have positive mass, temperature above 350 K, exact wall position
and component inventory. Deposition consumes both parcels; splash retains eight
children with IDs 3–10. Altered wall temperature rejects restart context.

These are application integration and restart checks. They do not establish
physical spray/film validity, full configuration coverage, selected-mechanism
qualification, installed-binary readiness or upstream solver parity. Rebalance
is requested, but actual ownership changes are not asserted by this checker.

## Completed results

Final CTest **1/1 passed**, 101.77 s total. Both deposition and splash complete
continuous/stopped/resumed runs and context-mismatch rejection. Final normalized
composite inventory errors are 9.1072982488782372e-18 (deposition) and
8.4324245413582203e-18 (splash), identical for continuous and resumed executions.
Gas, particle, film and budget comparisons pass at the existing 3e-12 relative
bound. Provenance, builds, Python syntax, project contract and diff whitespace
passed. The initial runtime failure described above was fixed; no build failures
occurred. Temporary diagnostic prints were removed before this final run.

## Reproduction

Artifacts: `build/spray-sdc-impact-app-checks/`; first failure retained in
`ctest-initial.log`, final test in `ctest.log`. Cases are under
`build/spray-amr-reactive-debug/spray-sdc-impact-restart/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc_impact_restart$' --output-on-failure
python3 -m py_compile tests/application/check_mpi_eb_impact.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
