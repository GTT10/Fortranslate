# Recursive LES in MPI AMR spray — 2026-09-25

Implementation `6862316`, based on `851db08`.
The existing `subcycling` flag now selects recursive LES half steps on both
sides of the reaction/transport/hydro split, inside the particle half steps.
False retains common LES stepping. The application passes the same flag to
the LES diffusion selector. Hydro and molecular-transport bounds remain
conservative finest-grid estimates. Existing interval retries remain active.

The coupled API agrees the effective subcycling flag across ranks before any
particle or LES work. LES boundary budgets and correction maxima join the
existing coupled ledger. Gas, temperature, parcels, component masses, optional
films and IDs still publish only after the entire interval succeeds. The
existing serialized `amr_spray` namelist already binds subcycling to checkpoint
identity; no new checkpoint field or format was introduced.

## Completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpiexec`), Ninja Debug.
Fixed tests-enabled application/unit and selected H2O2 tests-disabled application
built successfully. Selected execution was not run.

- `spray_mpi_amr_les_coupled_2`: passed, 63.42 seconds. Existing Vreman common
  and recursive cases exercise reacting continuation, periodic composite
  conservation and physical gas/parcel boundary budgets at the existing
  3e-12 scaled tolerance. New rank-disagreeing subcycling rejects before LES
  dispatch and preserves saved gas, temperature, particles, masses and IDs.
- `spray_mpi_amr_dynamic_subcycle`: final pass, 25.18 seconds (total 25.22).
  Gaussian Dynamic LES with recursive stepping, scheduled regrid and injections
  completes a two-rank run and a two-to-four-rank checkpoint continuation.
  Gas and parcel outputs agree at the existing rel_tol=2e-12, abs_tol=1e-25.
  Changing subcycling on restart rejects without final output. Existing changed
  filter/injection and malformed/incomplete-checkpoint checks also pass.

The initial application test failed after 15.65 seconds: its generated input
added `subcycling=.true.` before an existing false assignment, so both runs
actually used common stepping and the expected mode-change rejection did not
occur. The fixture now asserts exactly one false assignment and replaces it.
The corrected split input contains exactly one true assignment. The application
alone was rerun; the already-passed coupled test was not repeated. No numerical
tolerance was changed. Final evidence covers two tests; no skips.

Mechanism provenance, Python syntax, project contract and whitespace passed.
Compiler warnings include existing real-equality and injection-loop diagnostics.
MPI logs also contain host X authorization chatter; application return statuses
and assertions passed on the final run. This is focused coupled/application
integration evidence, not external solver parity, physical validation or a
full configuration matrix.

## Commands and artifacts

Logs: `build/les-mpi-spray-subcycle-checks/` (fixed/selected build logs,
initial `ctest.log`, final `ctest-application-final.log`).
Final application artifacts:
`build/spray-amr-reactive-debug/tests/spray_mpi_amr_dynamic_subcycle/restart_46_uvkx7/`.
Launcher environment overrides below are host-local workarounds.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_amr_spray_3d test_spray_mpi_amr_reactive --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_(dynamic_subcycle|les_coupled_2)$' --output-on-failure
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_dynamic_subcycle$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py tests/application/check_mpi_amr_spray_restart.py
python3 tools/check_project_contract.py --source .
git diff --check
```
