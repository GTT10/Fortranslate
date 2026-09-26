# Timed auxiliary sources for recursive EB SDC — 2026-09-27

Implementation `5a23dc7a`, based on `e3e913f0`, local `main`, version 0.246.0.
PeleC pin `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
Provenance preflight H2O2 source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
The focused unit fixture uses the existing elementary H2O2 mechanism.

## Interface and transaction

`advance_reactive_eb_tree_sdc_subcycle_3d` and
`advance_mpi_eb_tree_sdc_subcycle_3d` now accept paired optional
`auxiliary_source` and `auxiliary_change`. New public abstract interfaces are
`eb_sdc_recursive_source` and `mpi_eb_sdc_recursive_source`.

The serial callback is `(state,interval,stage,depth,time,rate,ok)`;
the MPI callback adds temperature after state and communicator before ok.
`time` is the level interval's start offset relative to the current outer call,
not absolute simulation time. It is unchanged across the old stage 0 and
new stages 1..iterations. Rate is per fluid volume per unit time. The existing
selected-level driver validates ownership/shapes, masks inactive levels and
solids, and computes the accepted old/final-new trapezoidal inventory.
The recursive driver accumulates those inventories separately from external
forcing and publishes only after all descendants and synchronization succeed.
Callbacks must keep any physical side effects in caller-private trial state.
The gas driver cannot undo arbitrary external callback mutations.

MPI checks callback/ledger pairing and rank presence before entering recursion.
The no-callback path explicitly omits both arguments. This is preparation for
recursive particle/film integration; it does not implement particle migration,
cross-level source routing or film scheduling within recursive intervals.

## Focused test scope

GNU Fortran 13.3.0, Ninja Debug, fixed-mechanism tests-enabled build
`build/spray-amr-reactive-debug`, OpenMPI 4.1.6 via `/usr/bin/mpifort` and
`/usr/bin/mpiexec`. Existing real-comparison and executable-stack warnings remain.

The new three-level/five-patch unit case uses two ranks, chemistry, molecular
transport and forcing, both distributed ownership and all patches on rank zero.
It checks 21 callbacks: three stages for each of seven level intervals, with
level sequence `[0,1,2,2,1,2,2]`, starts `[0,0,0,1,2,2,3]*dt/4` and intervals
`dt/2**depth`. The source depends on patch/cell, stage and interval start time.
Analytic integrated leaf source is `dt*(3-1/2**depth)` times its stage-zero
profile. The comparison detects intermediate-stage accumulation, duplicated
parent volume and wrong interval time/length.

Serial/MPI gas, temperature and reaction histories are compared exactly;
auxiliary inventories use 3e-14 relative and combined closure 3e-14 times
max(1,abs(initial inventory)). A callback rejection on the last descendant's
final stage, including on an idle rank, must retain gas, temperature, reaction
history, all ledgers, limiter and level-step counts. A no-source recursive
regression checks the existing path. Redistribution/source combinations are
not independently tested in this increment.

The initial build rejected an overlong Fortran call line, fixed by continuation.
The initial no-source run exposed a runtime check rejecting a disassociated
procedure pointer; explicit callback/no-callback calls fixed that path. The
initial auxiliary case passed. No tolerances or golden values were loosened.
No full matrix, installed audit, selected-mechanism rerun, temporal-order study,
upstream solver parity or physical validation was performed.

## Completed results

Final CTest **2/2 passed**, 116.59 s total: no-source recursive test 43.84 s,
auxiliary recursive test 72.70 s. Both ownership layouts produce exact
serial/MPI state, temperature and reaction history. Normalized composite
residual is **4.2613413117553738e-15** with auxiliary sources, and
3.9314133175908747e-15 without them. Analytic source inventory and late
collective rollback assertions pass. Provenance, build, Python syntax,
project contract and diff whitespace checks passed after the fixes above.
Final logs are `ctest-final.log` and `detail-final.log`.

## Reproduction

Artifacts: `build/sdc-recursive-source-checks/`; initial failures are retained
in `build.log`, `ctest-initial.log`, and `detail-initial.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_full_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_sdc_(aux_recursive|subcycle)_3d_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
