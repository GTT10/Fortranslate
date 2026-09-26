# Owner-local MPI tree StateRedist and common-step SDC — 2026-09-27

Implementation `8417dac2`, based on `45a09e1d` on local `main`.
PeleC pin: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`;
PelePhysics pin: `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`.
H2O2 source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Implemented scope

`redistribute_mpi_eb_tree_state_3d` uses replicated geometry and scalar
neighborhood metadata with owner-local signed component fields. Outward receiver
selection, overlap weights and connected-neighborhood fallback follow the
[serial composite operator](eb-tree-state-redistribution-20260927.md).
One sparse fetch supplies values for owned neighborhood means; another supplies
means contributing to owned recipients. Requests are deduplicated by leaf ID.
Neither gas fields nor neighborhood means are gathered into a full replicated
tree. Hidden parents retain their input, solid outputs are zero, and solid NaNs
are ignored. Boundary labels, target fraction and geometry agree collectively;
output is published only after successful finite calculations on all ranks.

The SDC adapter computes `(S(U+dt*R)-U)/dt`. Common-step MPI SDC applies it to
hydro and old/new molecular transport, matching serial. Reflux increments use
the linear state operator before EOS recovery and final restriction. Application
input now accepts `redistribution='state_redist'` and passes
`state_redist_target_volume_fraction` (default 0.5). Both remain part of restart
identity. Selected-level and recursive StateRedist remain rejected.

## Completed focused checks

GNU Fortran 13.3.0, Ninja Debug, fixed-mechanism tests-enabled build
`build/spray-amr-reactive-debug`; `/usr/bin/mpifort` and `/usr/bin/mpiexec`,
OpenMPI/OpenRTE 4.1.6. Existing compiler warnings and nonfatal launcher display
authorization diagnostics remain in logs.

- Provenance preflight passed. Fixed serial/MPI applications and MPI operator/
  coupled test targets built successfully.
- CTest **3/3 passed**, 18.51 s total: serial state operator 0.02 s,
  MPI coupled StateRedist SDC 17.90 s, MPI state operator 0.55 s.
- Root-only, left-refined, fully peer-refined and right-refined operator layouts
  cover outward receivers and connected fallback at coarse/fine interfaces.
  Both round-robin ownership and all patches on rank zero produce exactly the
  serial components, including constants. Composite signed-volume residuals
  range from 2.2204460492503131e-16 to 1.3322676295501878e-15.
- Rank-dependent target/boundaries and a nonfinite fluid value on one owner
  reject collectively and retain previous output. An idle rank participates.
- Three levels/five patches, two SDC iterations, reaction, molecular transport,
  hydro and external forcing match serial U, T and reaction estimates exactly
  for both ownership layouts. Inventories agree within 2e-12; limiter is exact.
  Composite relative residual is 4.0960116351235575e-15 for both layouts.
  Invalid viscosity after old transport and rank-dependent forcing preserve
  U/T/reaction estimates, inventories and limiter.
- Four application runs: serial continuous, MPI two-rank continuous, MPI
  two-rank stop after one step, and one-rank restart. The supplied StateRedist
  fixture advances to 2e-10 with dt=1e-10. Serial/MPI state CSVs are byte-identical;
  MPI continuous/restart state and budget CSVs are byte-identical. Maximum
  normalized budget residual is 2.975153346750058e-16.
- Changed-target restart and unsupported recursive input reject as expected.
- Python syntax, project contract and `git diff --check` passed.

An initial CMake registration duplicate and a subsequent Fortran line-length
error were corrected before the final successful build. Their logs are retained;
no tests were run from either failed build. Final CTest had no failures or skips.

Artifacts: `build/mpi-eb-tree-state-redist-checks/`, containing `provenance.log`,
`build.log`, `build-config-failure.log`, `build-line-length-failure.log`,
`ctest.log`, `detail.log`, `contract.log`, `run.py`, `app.log`, `summary.json`,
and isolated application inputs/logs/CSVs/checkpoint.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug \
  --target test_mpi_eb_tree_full_3d test_mpi_eb_tree_state_redistribution_3d \
  pelef_mpi_reactive_eb_tree_3d pelef_reactive_eb_tree_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 \
  OMPI_MCA_hwloc_base_binding_policy=none \
  ctest --test-dir build/spray-amr-reactive-debug \
  -R '^unit_(mpi_eb_tree_sdc_state_redist_3d_2|(mpi_)?eb_tree_state_redistribution_3d(_2)?)$' \
  --output-on-failure
python3 build/mpi-eb-tree-state-redist-checks/run.py
python3 -m py_compile tools/check_project_contract.py build/mpi-eb-tree-state-redist-checks/run.py
python3 tools/check_project_contract.py --source .
git diff --check
```

## Limits

No selected-mechanism build, full matrix, installed-binary audit, scaling study,
redistribution/regrid trajectory, temporal-order study, periodic StateRedist
case or upstream AMR scheduling comparison was run. Replicated geometry and
neighborhood searches remain a scaling limitation. Connected fallback targets
raw physical volume, not overlap-adjusted effective volume; these tiny-step
checks establish neither a general stable CFL nor physical flow validation.
MPI EB spray SDC and recursive StateRedist are still open.
