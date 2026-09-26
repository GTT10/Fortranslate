# Recursive spray SDC acceleration — 2026-09-27

Implementation `7a97fac36fd518ab37aa026ee26fee21cc1f99f8`, based on `94bf30220c4db0c85b3bd3471acaaa9064ea4cf1`,
local main, version 0.246.0.

## Change

Recursive MPI EB spray SDC now accepts uniform particle acceleration. The
application forwards its existing `acceleration` input to both gas and particles.
Each level applies a particle half kick before old-time exchange/drift and a
second half kick after final-iteration exchange. Intermediate SDC iterations
retain the cached old source and do not repeat acceleration.

Physical, virtual and ghost particles all receive acceleration. Only physical
particles contribute external momentum and kinetic work, using each half kick's
current mass and multiplicity. This ledger is separate from particle/gas exchange
and added to `source_change` only after the complete recursive transaction
succeeds. Private stage populations preserve rollback after later copy failure.
Acceleration must be finite, bounded by 1e6 m/s² per component, equal across ranks,
and unchanged across the stages of a level interval.

## Independently completed checks

Artifacts: `build/spray-recursive-acceleration-checks/`.
Existing fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`,
GNU Fortran 13.3.0, `/usr/bin/mpifort` and `/usr/bin/mpiexec` OpenMPI 4.1.6.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d test_spray_mpi_eb_sample test_spray_mpi_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_mpi_eb_sample_2|spray_mpi_eb_sdc_recursive_app)$' --output-on-failure
# After fixing the negative-test input override, rerun only that application test:
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc_recursive_app$' --output-on-failure
python3 -m py_compile tests/application/check_mpi_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Both focused tests passed on their final runs: population-stage test **1/1,
0.72 s**, and application test **1/1, 108.22 s** (108.27 s total for that rerun).
No skips. Provenance, build, Python syntax, project contract and whitespace checks
passed. Existing executable-stack linker warnings remain. Final application
artifacts are under
`build/spray-amr-reactive-debug/spray-sdc-recursive-app/eb_spray_d_s0lyn3/`.
The maximum scaled inventory residual was **2.126655773353596e-15** on both
rank layouts (`inventory.txt`); the 3e-12 criterion is unchanged.

The population-stage test uses acceleration `[1000,-2000,3000]` m/s² and checks
physical-only old/new force accounting, unchanged external work on intermediate
iterations, copy velocities after the final kick, changed rank-local acceleration
rejection, and late-copy failure followed by successful retry.

The application test runs three 10 ns outer intervals on 2 and 4 ranks with
chemistry, molecular transport, gas force density `[100,-200,300]` and acceleration
`[1000,-2000,3000]`. It compares rank layouts and checks the existing composite
inventory criterion of 3e-12, including external work. Negative launches cover
out-of-range acceleration, unsupported films and missing restart input.

The first combined test run passed the population test (0.72 s), then failed the
application checker (152.13 s). Its negative acceleration input had been inserted
before the valid fixture value, so the latter overrode it and the solver ran
successfully. The checker now appends overrides before the namelist terminator.
This was a test-input failure, not a numerical conservation failure. The first
run is retained in `ctest-input-override-failed.log` and
`LastTest-input-override-failed.log`; no tolerance was changed.

Recursive breakup and films remain unsupported. This check does not establish
accelerated rank-changing restart, selected-mechanism behavior, three-level
coupled accuracy, exact upstream parity or physical spray validation.
