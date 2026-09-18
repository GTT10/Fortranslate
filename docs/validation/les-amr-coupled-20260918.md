# MPI AMR reacting-spray LES coupling — 2026-09-18

Implementation `0a09c3616c8c5727328c6a08db74cf300b9de6ab`, based on
`ef7ba5d294fc5d88d4f362696576851702c82a72`.

`advance_mpi_amr_spray_reactive` accepts static LES settings and model. Two
common LES half steps surround the existing reaction/transport/hydro split,
inside the particle half steps. Smagorinsky, WALE and Vreman dispatch to the
existing distributed tree integrator. Defaults disable LES. Coefficients and
model agree collectively before any candidate work. LES boundary changes join
the gas/particle/film ledger; the reported maximum correction includes LES.
Private gas, temperature, parcel, component-mass and optional film candidates
publish only after the whole coupled interval succeeds.

Recursive gas stepping remains supported; LES half steps use a common interval.
This change does not wire application inputs, checkpoint identity or timestep
selection, and does not implement LES subcycling, dynamic LES or EB coupling.

## Completed verification

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; existing configurations unchanged.
Fixed tests-enabled `test_spray_mpi_amr_reactive` and selected H2O2 tests-disabled
`pelef_selected_spray_tree_runtime` built. Selected runtime was not run.

Initial CTest selection passed **3/3** in 172.20 seconds: LES coupling at
2 ranks (58.56 s), 4 ranks (55.93 s), and the default-disabled existing
`spray_mpi_amr_reactive_2` (57.66 s). The rank-disagreement fixture was then
strengthened to use a valid gas solver, rebuilt, and the two LES tests rerun:
**2/2 passed** in 115.33 seconds (59.22 s and 56.06 s). These are three distinct
tests, not five. Zero failed builds/tests or skipped tests occurred.

The coupled fixture uses Vreman 0.17, PrT/ScT 0.7 and CI 0.1 on a nested
three-node tree with owner-local gas and mixed water/peroxide parcels. It checks
common and recursive gas coupling, reacting continuation, periodic composite
conservation and slip-wall/outflow gas-plus-particle boundary accounting at the
existing `3e-12` relative tolerance. A later invalid gas solver rejects after
the first particle/LES halves while preserving caller fields, parcel IDs and
diagnostic sentinels. Rank-dependent CI is tested separately with a valid gas
solver. LES-on/off solution differences and film-plus-LES behavior were not
independently qualified in this fixture; other static closures were tested at
the underlying tree API in prior evidence.

Provenance, Python syntax, project contract and whitespace checks passed.
No tolerances changed. No full matrix, external solver parity, scaling or
physical validation was performed.

## Commands and artifacts

Logs: `build/les-amr-coupled-checks/`, containing `provenance.log`, `build.log`,
`build-selected.log`, `build-final.log`, `ctest.log` and `ctest-final.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_reactive --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_(reactive_2|les_coupled_[24])$' --output-on-failure
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_les_coupled_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
