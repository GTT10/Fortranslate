# Spray application LES selection — 2026-09-18

## Source and scope

Implementation: `afefec0e3241a8a757c61c279c39a058594edcbc`, based on
`0f425d6`. The regular serial and MPI spray drivers accept `les_model` with
Smagorinsky (default), WALE and Vreman choices. Fixed and selected frontends
share the drivers. Both SGS half steps receive the model, retaining existing
transactional diffusion rejection and step halving. Nondefault model names
enter checkpoint context; default context construction is unchanged. Invalid
models reject even with zero coefficient.

This is application plumbing and bounded restart evidence, not external PeleC
parity, turbulent-flow validation, nozzle validation or AMR/EB LES support.
The MPI fixture starts with uniform gas and develops particle-coupled gradients;
nontrivial prescribed-gradient closure checks are recorded separately in
[the MPI LES record](les-mpi-models-20260918.md). Selected executables were built,
but selected application runtime was not rerun in this increment.

## Independently completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug. Fixed tests-enabled build:
`build/spray-amr-reactive-debug`; selected H2O2 tests-disabled build:
`build/spray-amr-reactive-selected-debug`. No compiler or backend change.

- Mechanism provenance preflight passed.
- Fixed `pelef_spray_3d`, `pelef_mpi_spray_3d` and selected
  `pelef_spray_3d_selected`, `pelef_mpi_spray_3d_selected` built successfully.
- Serial existing Smagorinsky and new WALE/Vreman physical-boundary restart
  checks: **3/3 passed**, 6.63, 7.25 and 7.28 seconds. Wall/outflow cases retain
  exact gas/parcel CSV and final history replay, conservation gates and boundary
  behavior. Alternative model changes/removal and invalid names reject before
  history output.
- MPI WALE/Vreman checks: final **2/2 passed**, 11.54 seconds total. Two-rank
  continuous versus resumed gas/parcel records match exactly; 2-to-4-rank restart
  retains existing relative `2e-11`, absolute `1e-25` gates. Mass/energy budgets,
  overwrite refusal, model/context mismatch and invalid model rejection pass.
- Python syntax, project contract and `git diff --check` passed. No skipped tests.

Initial build failed on a new Fortran line exceeding 132 columns; wrapping it
fixed compilation. The initial five-test run passed all three serial tests but
failed both MPI checkers after successful restart: the checker expected the
word `checkpoint`, whereas the driver reports `Restart rejected: invalid data
or incompatible model/input context`. The assertion was corrected to that exact
existing diagnostic, and only the two MPI tests were rerun. No numerical gates
or application behavior were changed to address these test failures.

## Commands and artifacts

Logs: `build/les-application-model-checks/` (`provenance.log`, initial `build.log`,
`build-final.log`, `build-selected.log`, initial `ctest.log`,
`ctest-mpi-final.log`). Fixture inputs, outputs and per-invocation logs are under
`build/spray-amr-reactive-debug/tests/spray_{sgs,wale,vreman}_boundaries/` and
`tests/spray_mpi_{wale,vreman}_restart/` within the same build directory.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_spray_3d pelef_mpi_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_spray_3d_selected pelef_mpi_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_(sgs_boundary_restart|(wale|vreman)_boundary_restart|mpi_(wale|vreman)_restart)$' --output-on-failure
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_(wale|vreman)_restart$' --output-on-failure
python3 -m py_compile tools/check_spray_boundaries.py tests/spray/test_mpi_spray_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
