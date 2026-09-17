# Atomic EB gas/particle redistribution — 2026-09-18

Implementation: `0f00a99edfb04ad6ee39042cd30f90d90792008d`, with inventory-order fix
`c8b0cd3f492ca0b72d8f11df6a63e9270954b1d9`; original base `b8d0522`.

Optional `nodal_phi` extends `rebalance_mpi_amr_spray` to EB ownership. Mode
presence agrees across ranks; EB mode rejects periodic coordinates and film
arguments. Gas redistribution and EB particle/component migration operate on
private copies before publication. Rejected particle migration leaves the
original owner map, gas, temperature and particle population intact. EB geometry
validation still occurs if the proposed owner map is unchanged. Regular/film
paths retain their default behavior. EB composite integration now first reduces
one inventory column per patch (exactly one owner contributes), then sums those
columns in fixed tree order. This removes owner-dependent grouping of additions
in chemistry before/after inventories; it does not loosen a budget tolerance.

The fixed/selected MPI EB spray applications now connect the existing measured
gas work counters and balance controls to this joint publication path. Geometry
and hierarchy remain fixed. Counters measure gas kernels, not particle work;
this implementation is not a claim of optimal total-cost load balance or speedup.

## Completed focused checks

GNU Fortran 13.3.0 `/usr/bin/gfortran`, Open MPI 4.1.6 `/usr/bin/mpifort` and
`/usr/bin/mpiexec`, Debug/Ninja. Existing fixed tests-enabled and selected H2/O2
tests-disabled configurations were reused unchanged.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_balance test_spray_mpi_amr_balance pelef_mpi_eb_spray_3d test_mpi_eb_tree_full_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_mpi_eb_balance_[24]|spray_mpi_amr_balance_2|spray_mpi_eb_rebalance_application)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_mpi_eb_rebalance_application|unit_mpi_eb_tree_full_3d_2)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_mpi_eb_spray.py --exe build/spray-amr-reactive-selected-debug/pelef_mpi_reactive_eb_tree_3d_selected_spray --mpiexec /usr/bin/mpiexec --source . --work build/spray-mpi-eb-balance-checks/selected-integrals --selected --scheduled --restart --rebalance
python3 -m py_compile tests/application/check_mpi_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The library tests passed at two ranks (**0.47 s**) and four ranks (**0.49 s**),
plus the existing regular AMR balance regression (**0.53 s**). After the inventory-order fix, the existing MPI EB full gas test passed in
**42.60 s**, and the final fixed application test passed in **54.99 s** (the two final tests total **97.62 s**). Selected
application execution passed.
Provenance, Python syntax, contract and whitespace passed; no skips.

After real redistribution was exercised, the fixed continuation test exposed a
chemical-ledger difference of about 2.3e-24 (ledger about -1.706e-15), caused by
rank-dependent summation grouping. The fixed-order EB inventory reduction above
addresses this numerical reproducibility issue. No tolerance was changed.

Two initial application fixture variants failed the assertion that a move
occurred: with two patches and at least two ranks, swapping owners cannot improve
the peak load. Solver runs themselves succeeded. The final fixture uses a root
and four nonoverlapping ratio-two xy quadrants ordered low-y/high-y/low-y/high-y,
with different fluid work near the y=1.25 mm EB. Both fixed and selected runs
confirm actual `Rebalanced at step:` messages. This fixes the fixture rather
than requiring an already balanced two-patch simulation to move.

The library fixture starts both patches on rank N-1 with costs [1,100]. It first
rejects a zero particle capacity after private gas redistribution and checks
owner, gas, temperature and component rollback. With adequate capacity it moves
the fine patch and its mixed parcel to rank zero, preserving active gas values,
covered zero states/temperature sentinels and particle/component identity.
The application fixture compares all gas/particle CSV values at two/four ranks,
then checks two-rank checkpoint to four-rank continuation against uninterrupted
four-rank output, including pending injection events and ledgers. Existing
context/truncation and immutable-publication checks also run.

Artifacts: `build/spray-mpi-eb-balance-checks/` contains initial `build-final.log`,
`build-selected.log`, `ctest.log`, intermediate `ctest-application-final.log`, and
final `build-integrals.log`, `build-selected-integrals.log`, `ctest-integrals.log`,
`selected-integrals.log`.
The first CTest log retains the three successful library/regression results and
the initial application failure. Final selected outputs are
`selected-integrals/eb_spray_jk_n13ys/`; the CTest LastTest log records the final
fixed fixture directory. Earlier failed fixtures are retained for diagnosis.
No complete matrix, performance/scaling qualification, regridding, external
solver parity or physical validation was run in this increment.
