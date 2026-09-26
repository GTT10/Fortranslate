# Source-level EB redistribution transfers — 2026-09-27

Implementation `fdc880e4`, based on `77ea15d8` on local `main`.
PeleC pin: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

The serial and MPI libraries now expose
`redistribute_sdc_eb_tree_level_rate_3d` and
`redistribute_sdc_mpi_eb_tree_level_rate_3d`. They decompose the linear composite
flux-redistribution operator by the level supplying a rate. Rates on all other
levels are zeroed before redistribution; their payloads are not used.

The two outputs are:

- `rate`: redistributed rate on the selected source level, zero elsewhere.
- `transfer`: dt times the redistributed rate on other levels, zero on the
  selected level. This is an increment per fluid volume for a caller-owned
  interval register, not a cell-integrated mass.

Thus the fluid-volume integral of `dt*rate + transfer` equals the selected
input's integrated contribution. Summing `rate + transfer/dt` over source
levels reconstructs the original composite rate operator. Both outputs publish
only after validation and calculation succeed. The MPI form retains owner-local
storage and collectively agrees source depth and positive finite interval before
exchange. Other-level NaNs are ignored; selected fluid values retain the existing
finite-value requirement. Transfer multiplication is guarded against overflow.

These routines do not advance gas, accumulate a persistent interval register or
change application dispatch. Recursive SDC redistribution remains disabled until
its stage transfers and synchronization are connected.

## Completed focused checks

GNU Fortran 13.3.0, `/usr/bin/mpifort`, `/usr/bin/mpiexec`
(OpenMPI/OpenRTE 4.1.6), Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`; existing warnings remain.

- Mechanism provenance and both touched test targets built successfully.
- CTest: 2/2 passed, 0.63 s total. Serial operator test 0.03 s; two-rank MPI
  operator test 0.55 s.
- Root-only, partially refined and peer-refined layouts exercise each source
  level, with periodic and coarse/fine impulse cases. Intervals are
  0.25/(2**depth). Each level's integrated contribution closes within 3e-14
  absolute discrepancy; reconstruction agrees within 2e-14 per component.
- Other-level inputs are deliberately NaN. Same-level deferred outputs and
  other-level immediate outputs are exactly zero.
- MPI rates and transfers match serial exactly for round-robin ownership and
  all patches on rank zero, retaining owner-local allocation.
- Negative serial intervals, rank-dependent depths and rank-dependent MPI
  intervals reject without changing either output.
- Python syntax, project contract and `git diff --check` passed. No failures,
  skipped tests or interrupted runs were counted.

Artifacts: `build/sdc-eb-redist-level-checks/`, including `provenance.log`,
`build.log`, `ctest.log`, `detail.log` and `contract.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug \
  --target test_eb_tree_redistribution_3d test_mpi_eb_tree_redistribution_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 \
  OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug \
  -R '^unit_(mpi_)?eb_tree_redistribution_3d(_2)?$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

This establishes a spatial decomposition and transfer accounting, not recursive
time integration, temporal order, a larger stable step or upstream AMR
redistribution parity. No application trajectory, selected build, installed
binary, full matrix or physical validation was run in this increment.
