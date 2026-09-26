# SDC reaction history during spray load balancing — 2026-09-27

Implementation `1127c62dbe4e477da0ddd0dbe5ad64a273715f12`, based on `18af1832`,
local main, version 0.246.0.

## Fix

The MPI EB gas-only application already moved its SDC reaction-source history
when gas ownership changed. Its spray branch instead called the gas/particle/
film rebalance helper without the history, leaving that owner-local field on
old ranks. The next SDC interval or checkpoint could therefore encounter
inconsistent ownership.

`rebalance_mpi_amr_spray` now accepts optional `reaction_source`, with nvar-4
components. Presence is agreed across ranks before the branch. It uses the
existing gas-plus-field rebalance helper on private gas, temperature, history
and owner arrays. The history publishes only after physical particle and any
film migration have succeeded. The application passes its optional SDC history
and updates the checkpoint record after a successful transaction. Both common
and recursive SDC use this application branch; non-SDC callers omit the field.

## Independently completed checks

Artifacts: `build/spray-sdc-rebalance-history-checks/`.
Existing fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`,
GNU Fortran 13.3.0, `/usr/bin/mpifort` and `/usr/bin/mpiexec` OpenMPI 4.1.6.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_balance pelef_mpi_eb_spray_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_balance_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final tests **2/2 passed**, 0.50 s on 2 ranks and 0.51 s on 4 ranks, 1.06 s
total. No test failures or skips. Provenance, final build (including application
wiring), Python syntax, contract and whitespace checks passed.

The existing two-level planar EB test now assigns a distinct exact value to
every reaction-history component/cell, including covered cells. It checks:

- reaction history is allocated only on the accepted new owner;
- every transferred history value is identical to its original value;
- missing history on one rank rejects collectively;
- particle-capacity and later film-capacity rejection retain gas, temperature,
  physical parcels, films, owners and reaction history.

The initial build failed because the extended declaration exceeded Fortran's
line-length limit. The declaration was wrapped, then both final builds passed;
`build.log`, `build-final.log`, and `build-test.log` retain the sequence.
`ctest.log` and `LastTest.log` contain the final tests. Existing build warnings
remain. No end-to-end timed application rebalance, post-rebalance restart,
selected-mechanism matrix or physical validation was run in this check.
