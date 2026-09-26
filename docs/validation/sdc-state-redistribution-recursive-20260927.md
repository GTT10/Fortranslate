# Conservative StateRedist transfers and recursive SDC — 2026-09-27

Implementation `acd7f9dc`, based on `c084efbf` on local `main`.
PeleC pin: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`;
PelePhysics pin: `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`;
AMReX pin: `1546c1ccbf41c704a33790f4bf445a2fcee1372e`.
H2O2 source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## State transfer decomposition

StateRedist changes the underlying state even at zero input rate. Simply
selecting a level's rate and applying the linear redistribution misses this
baseline displacement. The serial/MPI composite operators now accept optional
`correction_depth`: return signed displacements from neighborhoods rooted at
that depth, with zero hidden-parent/solid entries, instead of a full state.
Neighborhoods and weights are unchanged from the common-step operator.

Let `S` be composite linear state redistribution and `D_l(U)` the sum of
neighborhood displacements rooted on level `l`. For each neighborhood, its
member's contribution is its partition coefficient times `(mean - U_member)`.
Every neighborhood's fluid-volume integral is zero. Thus every level's baseline
correction is conservative and `sum_l D_l(U) = S(U)-U`, up to roundoff.
Outward receivers and connected-neighborhood fallback both use this accounting.

The new serial/MPI selected-level adapters return
`S(P_l R) + D_l(U)/dt` on the selected level, where `P_l` selects input rates on
level `l`. Other levels receive separate integrated increments
`dt*S(P_l R) + D_l(U)`. Summing all levels' rates and transfers at a common state
reconstructs `(S(U+dt*R)-U)/dt`. Zeroing unselected rates does not remove the
baseline state used in neighborhood means. MPI exchanges remain sparse and
owner-local; correction-depth selection is collectively validated.

Selected-level SDC now accepts `state_redist`, requires transfer output, and
publishes the final hydro plus half old/new transport transfers. Intermediate
SDC iterations replace these transfers. Recursive serial/MPI drivers use the
existing interval registers: a transfer from level `s` to level `r` waits in
bucket `min(s,r)` until that interval closes. Reflux increments pass through the
linear state operator before accepted register corrections, EOS recovery and
restriction. State, temperature, reaction estimates and ledgers stay private
until the entire root interval succeeds.

Applications now allow `redistribution='state_redist',subcycling=.true.` and pass
the target fraction through recursive dispatch. Method, target and time-step
configuration remain part of checkpoint identity.

## Reference boundary

The pinned `Src/EB/AMReX_EB_StateRedistribute.cpp` was inspected using `git show`
in `/mnt/sn850/Projects/gtt2/external/PeleC/Submodules/PelePhysics/Submodules/amrex`.
It forms weighted neighborhood means and tracks individual coarse/fine movements
in `MLStateRedistribute`. This increment decomposes PeleF's existing composite
leaf operator into conservative neighborhood corrections. It does **not**
establish equivalence to AMReX's level-local ghost/covered-cell movement ledger
or its asynchronous scheduling. The captured reference file is in the artifacts.

## Completed focused checks

GNU Fortran 13.3.0, Ninja Debug, fixed-mechanism build
`build/spray-amr-reactive-debug`; `/usr/bin/mpifort` and `/usr/bin/mpiexec`,
OpenMPI/OpenRTE 4.1.6. Existing compiler warnings and nonfatal launcher display
authorization diagnostics remain in logs.

- Provenance passed. Serial/MPI applications, serial/MPI state operator tests
  and MPI coupled SDC test built successfully.
- CTest **3/3 passed**, 44.23 s total: serial state operator 0.03 s, MPI recursive
  StateRedist SDC 43.56 s, MPI state operator 0.58 s. No failures or skips.
- Four operator layouts cover root-only, left refinement, fully refined peers
  and right refinement, including coarse/fine outward receivers and connected
  fallback. Per-level baseline corrections conserve signed fluid-volume
  integrals within 3e-14 absolute. Their sum reconstructs the full state
  displacement within 3e-14; rates plus transfers reconstruct the full affine
  source operator within 3e-14. Each level's combined update has the inventory
  of its selected raw rate within 3e-14.
- MPI correction, selected rate and transfer fields equal serial exactly under
  round-robin and all-rank-zero ownership, including an idle rank. Unselected
  rates and selected-level transfers are zero. Existing invalid target/boundary
  and one-owner nonfinite-input rollback checks also pass.
- Three levels/five patches, two iterations, hydro, transport, reactions and
  external force: recursive counts `[1,2,4]`; exact serial/MPI U, T and reaction
  estimates for both ownership layouts; exact limiter and inventory agreement
  within 2e-12. Composite relative residual is 3.9314133175908747e-15 for both
  layouts. Invalid viscosity preserves U/T/reaction estimates, step counts,
  inventories and limiter after old transport has run.
- Four application runs use `sdc_state_redist.nml` with `subcycling=.true.`,
  dt=1e-10 and final time 2e-10: serial continuous, MPI two-rank continuous,
  MPI two-rank stop after one step, and MPI one-rank restart. Serial/MPI state
  CSVs are byte-identical; MPI continuous/restart state and budget CSVs are
  byte-identical. Maximum normalized budget residual is
  4.1016292491558864e-16. Changed-target and changed-method restart reject.
- Python syntax, project contract and `git diff --check` passed.

Artifacts: `build/sdc-state-redist-recursive-checks/`: `provenance.log`,
`build.log`, `ctest.log`, `detail.log`, `contract.log`, `toolchain.log`,
`reference-state-redistribute.cpp`, `run.py`, `app.log`, `summary.json`, and
isolated application inputs/logs/CSVs/checkpoint.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug \
  --target test_eb_tree_state_redistribution_3d test_mpi_eb_tree_state_redistribution_3d \
  test_mpi_eb_tree_full_3d pelef_reactive_eb_tree_3d pelef_mpi_reactive_eb_tree_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 \
  OMPI_MCA_hwloc_base_binding_policy=none \
  ctest --test-dir build/spray-amr-reactive-debug \
  -R '^unit_(mpi_eb_tree_sdc_state_recursive_3d_2|(mpi_)?eb_tree_state_redistribution_3d(_2)?)$' \
  --output-on-failure
python3 build/sdc-state-redist-recursive-checks/run.py
python3 -m py_compile tools/check_project_contract.py build/sdc-state-redist-recursive-checks/run.py
python3 tools/check_project_contract.py --source .
git diff --check
```

## Limits

No temporal-order, stable-CFL, upstream multilevel scheduling, physical flow,
redistribution/regrid trajectory, selected-mechanism build, installed-binary or
full-matrix qualification was run. Conservation and serial/MPI agreement in these
bounded cases do not replace those evidence classes. Neighborhood metadata remains
replicated. MPI EB spray SDC remains unconnected.
