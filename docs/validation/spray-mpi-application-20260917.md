# Fixed MPI spray executable — 2026-09-17

Parent `fd1ad3e` plus this change. Build/test/preflight logs and source hashes:
`build/spray-mpi-application-checks/`. Application artifacts remain in fresh
`run-*` directories under
`build/spray-mpi-migration-gnu-debug/tests/spray_mpi_application_{2,4}/`.

The executable reads a root-broadcast namelist, injects a two-component cone,
migrates particles, chooses a global hydro/transport time step, retries rejected
whole coupled steps, and writes inventory history and rank-local gas/particle CSVs.
It refuses existing outputs. This is fixed-mechanism, periodic x-slab operation;
CSV output is not a restart checkpoint.

GNU Fortran 13.3.0 Debug, Open MPI 4.1.6, unchanged existing MPI configuration
from `spray-mpi-migration-20260916.md`:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target pelef_mpi_spray_3d --parallel 8
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_application_' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
python3 -m py_compile tests/spray/test_mpi_spray_application.py
git diff --check
```

2/2 application entries pass (2 and 4 ranks), 2.62 s, no failures/skips. Each
runs two coupled steps to 2e-8 s, checks relative mass/energy drift below 2e-11,
all 128 global cells exactly once, finite positive density/temperature, nonnegative
species and closure, four unique particle IDs, component inventories and particle
heating. A repeated invocation fails and preserves the history bytes. Provenance,
contract, Python syntax and whitespace checks pass. Runtime kernels retain Debug
bounds/FPE checks; the new program target uses the compiler's standard Debug flags.

No full matrix, selected-mechanism MPI executable, checkpoint/restart, scheduled
injection, physical validation or production scaling is claimed.
