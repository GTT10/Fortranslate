# MPI EB source coupling with KH children — 2026-09-18

Implementation `63e0d59`, base `ec059ec`.

The MPI EB mixture source now admits KH child mode with a shared ID watermark.
Before release, it checks the source-updated parent trajectory against the
original nodal surface. Crossing parents skip aerodynamic release for that
substep, preserving wall precedence. An optional eligibility mask in the
collective KH helper leaves excluded parent/pending state unchanged and returns
zero release energy for them. Unhit parents can create four children; the
transverse-energy debit joins gas feedback at the parent source location.
Parents and children then share specular or impact flight and final migration.
The whole source interval remains transactional.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpifort`, `/usr/bin/mpiexec`),
Debug/Ninja. Existing fixed tests-enabled and selected H2O2 tests-disabled builds
were reused; the EB impact source test retains runtime checks and floating traps.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_impact_step test_spray_mpi_amr_step --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_(eb_impact_step|amr_step)_[24]$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_impact_step --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_impact_step_[24]$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The initial run passed regular AMR **2/2**, while both new EB cases failed an
owner expectation: the parent endpoint lay exactly on the fine/coarse interface,
so placing every child on the coarse side was not guaranteed. Moving the fixture
start 5 micrometres toward that interface gives an unambiguous crossing. Final
EB CTest **2/2 passed**, 1.24 s total (0.59/0.60 s for 2/4 ranks). Thus four
distinct tests have passing results; the two earlier failures remain recorded.
No skips. Selected library build, provenance, Python syntax, contract and
whitespace passed. No selected execution was run.

The wall fixture enables KH children during a 200 m/s impact but still creates
only the four splash children and one film, with no preceding KH release.
Existing empirical-model rejection and whole-interval rollback checks remain.
The no-crossing fixture uses a mixed water/peroxide parent at 1000 m/s,
1e-8 s interval and reference diameter cubed 1e-25 m3. It creates four KH
children with watermark five, preserves component fractions, advances flight
and migrates to the exposed coarse owner. Both impact-enabled and specular
source calls conserve composite gas/particle/wall inventory within relative
3e-12, scale floor 1e-30. A capacity-one/substep-one attempt retains original
particle payload, watermark and ledger on failure.

Artifacts: `build/spray-mpi-eb-khrt-checks/{build.log,build-final.log,ctest.log,ctest-final.log,build-selected.log,provenance.log,contract.log}`.

This connects the source operator. MPI reactive specular shared-ID handling and
application KH mode remain to be connected. The operator remains source-first
breakup followed by flight, not globally time-ordered aerodynamic breakup.
No broad matrix, installed audit, external parity, scaling or physical breakup
validation was run.
