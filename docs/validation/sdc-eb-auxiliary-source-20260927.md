# Old/new auxiliary sources in EB tree SDC — 2026-09-27

Implementation `468d6d29`, based on `7585a13c` on local `main`.
PeleC pin: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`;
PelePhysics pin: `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`.
H2O2 source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Gas-side stage connection

Serial and MPI `advance_*eb_tree_sdc_3d` now accept paired optional arguments
`auxiliary_source` and `auxiliary_change`. The source procedure receives the old
state at stage 0, and the provisional new state at stages 1 through `iterations`.
Its inputs include the interval and selected depth (`-1` for a common step).
The MPI procedure additionally receives owner-local temperature and the driver
communicator. The driver supplies an allocated zero source field with matching
ownership; the procedure fills conserved rates per fluid volume per time.

The rates enter CTU prediction, old/new predictor export and the nonreacting
forcing passed to the SDC chemistry solve. Subsequent iterations replace the
new source estimate. The accepted auxiliary inventory is
`dt/2 * integral(old_rate + final_new_rate)`; intermediate source estimates are
not accumulated. `source_change` retains only the original prescribed-force
inventory. `auxiliary_change` separately accounts for mass, momentum, energy and
species supplied by the auxiliary family, which is necessary when gas/particle
exchange must cancel inside a combined-system budget.

The output procedure and inventory must be supplied together. MPI checks their
presence collectively before branching. Source shapes/ownership are validated;
solid and inactive-level source values are zeroed before finite checks. An
error returned on any rank rejects collectively. State, temperature, reaction
estimates, predictors and all inventories publish only after the whole driver
succeeds. A callback's own external side effects are **not** automatically rolled
back: a particle caller must keep trial parcels, films and IDs private until
success.

The hook is implemented in the common/selected-level driver. Recursive wrappers
do not yet expose it. This increment supplies the gas-side consumer for the
[independent particle source export](spray-sdc-source-export-20260927.md).
Particle kick/drift scheduling, particle rollback around the complete SDC driver
and the spray application's SDC selector remain to be connected.

## Reference inspection

Pinned `Source/Sources.cpp` in `/mnt/sn850/Projects/gtt2/external/PeleC` routes
old/new spray sources to `particleMKD`/`particleMK`; the particle routines were
inspected in the preceding source-export increment. Here those source-family
semantics motivate a distinct gas-side callback and inventory. No full particle
trajectory or upstream SDC parity claim is made.

## Completed focused checks

GNU Fortran 13.3.0, Ninja Debug, fixed-mechanism tests-enabled build
`build/spray-amr-reactive-debug`; `/usr/bin/mpifort` and `/usr/bin/mpiexec`,
OpenMPI/OpenRTE 4.1.6. Existing compiler/display diagnostics remain. Passing
internal Fortran test procedures also produces the linker's executable-stack
warning for the test executable.

- Mechanism provenance and `test_mpi_eb_tree_full_3d` build passed.
- `unit_mpi_eb_tree_sdc_aux_source_3d_2`: **1/1 passed**, 43.79 s test,
  43.84 s total; no failures or skips.
- Three levels/five patches, two iterations, dt=1e-10, molecular transport,
  hydro, reactions and prescribed force. A separate spatially varying source
  injects density, momentum, energy and one species with stage multipliers
  1, 2 and 3. Both serial and MPI observe exactly stages `[0,1,2]`.
- The auxiliary inventory agrees within 3e-14 relative with the independently
  integrated base profile times `2*dt`, proving use of old/final-new weights
  rather than accumulation of the intermediate stage.
- Round-robin and all-rank-zero ownership produce exact serial/MPI gas,
  temperature and reaction estimates, exact limiter, other inventories within
  2e-12 normalized discrepancy, and auxiliary inventories within 3e-14 relative.
  Composite normalized residual is **4.2012936932482202e-15** for both layouts,
  below the unchanged 3e-14 bound.
- Solid source NaNs are ignored. A callback error on one rank (including an
  idle rank), a fluid NaN on one owner, and a missing owner field all reject at
  the final source stage after two successful calls. Gas, temperature, reaction
  estimates, auxiliary inventory, existing inventories and limiter retain
  their input values.
- Python syntax, project contract and `git diff --check` passed.

Artifacts: `build/sdc-eb-aux-source-checks/`: `provenance.log`, `build.log`,
`ctest.log`, `detail.log`, `contract.log`, `toolchain.log` and
`reference-sources.cpp`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug \
  --target test_mpi_eb_tree_full_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 \
  OMPI_MCA_hwloc_base_binding_policy=none \
  ctest --test-dir build/spray-amr-reactive-debug \
  -R '^unit_mpi_eb_tree_sdc_aux_source_3d_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

## Limits

This is a prescribed-source common-step check, not an evolving-particle SDC
calculation. No source-enabled selected-level test, recursive source schedule,
application/restart, selected-mechanism build, source/redistribution combination,
full matrix, installed-binary audit, temporal-order study or physical validation
was run.
