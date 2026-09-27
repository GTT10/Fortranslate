# Film-bearing recursive MPI restart — 2026-09-27

Regression `e3fd680765cf4e28d51a047a9b834fb252231b4f`, based on `ec6348a2`,
local main, version 0.246.0. Solver source is unchanged by this increment.

## Check

The recursive impact application checker now has `--restart`. It runs a
four-rank continuous calculation for two root intervals, a two-rank calculation
stopped after the first interval, and a four-rank continuation from that archive.
Gas chemistry, molecular transport, particle exchange, splash and stationary
film evolution remain enabled throughout.

The stopped particle archive must contain a nonzero film count and an advanced
particle-ID watermark. Continuation must finish the second interval. Final gas,
sorted particle and sorted film CSV records must exactly match the continuous
calculation. All accepted states must have finite values, positive masses,
closed component inventories and unique particle IDs. Combined inventories use
the existing `3e-12` criterion, with the established mass/energy scale for
near-zero momentum. All stored inventory columns are also compared between
continuous and restarted output. Reading the checkpoint must leave the manifest,
gas file and particle/film file byte-identical.

## Independently completed checks

Artifacts: `build/spray-recursive-film-restart-checks/`. Case outputs:
`build/spray-amr-reactive-debug/spray-sdc-recursive-impact-restart/impact_25juorid/`.
Existing fixed-mechanism Ninja Debug build; GNU Fortran 13.3.0 and OpenMPI 4.1.6
(`/usr/bin/mpifort`, `/usr/bin/mpiexec`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tests/application/check_recursive_eb_spray_impacts.py tools/check_project_contract.py
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc_recursive_impact_restart$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 CTest passed**, comprising three application runs, 72.35 s test time and
72.39 s total. Provenance, build, Python syntax, project contract and whitespace
checks passed. No failures or skips. Logs are `build.log`, `ctest.log` and
preserved `LastTest.log`.

Continuous and restarted calculations both finish at `2e-8 s` with four
particles and one film. Gas, particle and film CSVs are **exactly equal**.
The maximum normalized inventory residual is **1.4082476556244309e-15** in
both outputs. Archive bytes are unchanged after continuation. No solver repair
was needed.

This is fixed-topology, fixed-mechanism film-bearing continuation over two root
intervals. Film-bearing regrid/rebalance, long-time evolution, selected
mechanisms, physical validation and full external-solver parity were not run.
