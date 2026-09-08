# Contributing to PeleF

## A small, testable change

Choose one issue and its acceptance conditions. Use one branch and one PR;
after integration start from an updated `main`, not an old implementation
branch. Preserve uncommitted user work. Commit and push only within the user's
authorized task. Do not merge failing checks or claim unrelated work complete.
The temporary exception is the existing 0.245.0 integration PR #102: finish its
checkpoint repairs before adding another feature stack.

## Local checks

From the repository root:

```sh
python3 -m unittest discover -s tests/python -p test_mechanism_provenance.py -v
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
python3 -m compileall -q tools tests/python
git diff --check
cmake -S . -B build/debug -DCMAKE_BUILD_TYPE=Debug
cmake --build build/debug --parallel 2
ctest --test-dir build/debug --output-on-failure
```

Use `ctest -R` for a focused iteration, then report the regex and count rather
than calling that run the full suite. A fresh Release directory is separate
from Debug. Select parallelism to fit memory and MPI oversubscription.

## Optional configurations

Install Cantera into a virtual environment rather than the system Python:

```sh
python3 -m venv build/cantera-env
build/cantera-env/bin/python -m pip install cantera==3.2.0
build/cantera-env/bin/python tools/check_mechanism_provenance.py \
  --bundle mechanisms/h2o2_full.json --regenerate
cmake -S . -B build/debug-cantera -DCMAKE_BUILD_TYPE=Debug \
  -DPELEF_ENABLE_CANTERA_REFERENCE=ON \
  -DPython3_EXECUTABLE="$PWD/build/cantera-env/bin/python"
cmake --build build/debug-cantera --parallel 2
ctest --test-dir build/debug-cantera --output-on-failure
```

MPI uses `-DPELEF_ENABLE_MPI=ON` in a new directory. Supply matching
`CMAKE_Fortran_COMPILER`, `MPI_Fortran_COMPILER`, and `MPIEXEC_EXECUTABLE` when
more than one toolchain is installed. Selected-mechanism builds add
`-DPELEF_MECHANISM_BUNDLE=mechanisms/h2o2_full.json` and
`-DPELEF_MECHANISM_SOURCE=mechanisms/h2o2_cantera.yaml`.

Optional CVODE adds `-DPELEF_ENABLE_SUNDIALS=ON` and a `CMAKE_PREFIX_PATH` to
SUNDIALS 7.2.0 with official static Fortran modules and double precision.
This is presently a selected 0D backend, not general stiff CFD chemistry.
Do not treat an unavailable dependency as a passed or intentionally disabled
qualification gate.

## Evidence and CI

Record the commit and tree (`git rev-parse HEAD` and `git rev-parse HEAD^{tree}`),
local dirty state, compiler/dependency versions, full commands and test results.
A commit name alone does not identify a dirty worktree. Use CTest JUnit output
and retain `Testing/Temporary/LastTest.log`, `CMakeCache.txt`, and any install
manifest alongside the report. Label failures, skips and interrupted runs.

The mechanism-integrity workflow compares exact YAML bytes and complete bundle
regeneration before numerical compilation. It never repairs inputs in CI.
The MPI installed audit is in `tools/ci_mpi_installed_smoke.sh`; run it only
with a fresh matching `build-mpi` / `install-mpi` pair. Its 42-file contract is
17 fixed serial + 10 selected serial + 11 fixed MPI + 4 selected MPI programs.
Tests-disabled release qualification is distinct from this tests-enabled CI
install and from CTest's private staged prefixes.

PRs must explain numerical impact, supported scope, test evidence, and remaining
limits. Retain historical evidence unchanged and add new dated results; do not
silently relabel an older source's successful run as the current source.
