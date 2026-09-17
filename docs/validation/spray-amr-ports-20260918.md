# Multiport MPI AMR spray injection — 2026-09-18

Implementation: `3488a4092ad29cba1e7661f55d60f4538478483c`, based on
`cab5acc`. The shared fixed/selected driver supports up to 16 ports with
separate center, axis, speed, cone angle, diameter, temperature, radius,
composition and existing nozzle profile/swirl/spread settings. Components use
a shared liquid table. Events select a port; equal-time events execute in array
order after any coincident regrid, and checkpoints follow all events at that
time. Each event retains private injection/migration and ID/ledger publication.
Newborn `material_id` records the port and is preserved through migration and
restart. Default port-count zero maps the previous scalar settings to port one.

Port inputs are included in the numerical context. Saved contexts from earlier
application builds do not acquire automatic compatibility with newly serialized
controls; the library checkpoint payload format remains V1.

## Completed focused verification

Reused GNU Fortran 13.3.0 Debug/Ninja fixed tests-enabled and selected H2/O2
tests-disabled builds, `/usr/bin/gfortran`, Open MPI `/usr/bin/mpifort` and
`/usr/bin/mpiexec`. Compiler/MPI/backend configuration did not change.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_amr_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_spray_3d_selected --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_(ports|application)$' --output-on-failure --output-log build/spray-amr-ports-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_mpi_amr_spray_restart.py --exe build/spray-amr-reactive-selected-debug/pelef_mpi_amr_spray_3d_selected --mpiexec /usr/bin/mpiexec --source /home/gtt13/projects/pelec --work /home/gtt13/projects/pelec/build/spray-amr-ports-checks/selected --selected --multiport
python3 -m py_compile tests/application/check_mpi_amr_spray_restart.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Fixed CTest **2/2 passed**, no skips: multiport restart **4.96 s**, existing
single-port application regression **4.28 s**, **9.28 s** total. Selected
multiport restart checker also passed.

The two-port fixture uses different locations, opposing axes, speeds, diameters,
temperatures and water/peroxide fractions, with finite solid/hollow outlets and
swirl. Four events occur at 0, 10, 10 and 15 ns. The checkpoint at 10 ns follows
coarsening and both coincident events. Two-rank save/four-rank continuation
matches uninterrupted gas/particle CSVs at relative 2e-12 / absolute 1e-25.
Checks confirm four consumed events, 16 surviving parcels, eight records per
port ID, separated outlet regions and the expected ordering of component masses
for each port. Changed port speed rejects context; truncated particles and
overwrite attempts also reject as expected. Single-port mass/component budgets
and rank parity retain their existing regression coverage.

No build/test attempt failed in this increment; the rejection subprocesses are
expected assertions. Provenance, Python syntax, project contract and whitespace
checks passed. Compiler warnings remain. No new independent multiport energy
budget, active-breakup restart, installed audit, convergence, scaling or physical
nozzle/spray validation was performed. Gas chemistry/transport are disabled in
these application fixtures.

Logs: `build/spray-amr-ports-checks/{build,build-final,build-selected,ctest,test-selected}.log`.
Fixed multiport artifacts: `build/spray-amr-reactive-debug/tests/spray_mpi_amr_ports/restart_*`.
Selected artifacts: `build/spray-amr-ports-checks/selected/restart_*`.
