# MPI breakup application and restart — 2026-09-17

Source: parent `72a2524` plus this change. Hashes and build/preflight/test logs in
`build/spray-mpi-breakup-app-checks/`. Existing GNU Fortran 13.3.0 Debug / Open MPI
4.1.6 configuration.

The MPI application now accepts representative ETAB/KH-RT breakup, common liquid
surface tension/viscosity and KH-RT constants. Input validity and checkpoint
identity include these settings. Particle CSVs expose deformation, deformation
rate and RT time. The MPI checkpoint API accepts optional breakup/RT payload
flags, checks collective agreement, and uses the existing serial payload schemas.
Migration already retains these fields. Defaults preserve disabled-breakup context
and schema. Child shedding and port-specific breakup properties remain open.

Completed commands:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target pelef_mpi_spray_3d test_spray_mpi --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_(application|migration)_[24]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
python3 -m py_compile tests/spray/test_mpi_spray_application.py
git diff --check
```

Four focused entries passed (4/4, 62.16 s; no failures/skips). The two application
entries now run ETAB/KH-RT with 1000 m/s injection through 2e-8 s, checkpoint at
1e-8 s and continue on changed ranks (2↔4). Final gas and particles, including
breakup state, match continuous output at relative 2e-11 / absolute 1e-25. ETAB
has positive deformation; KH-RT increases multiplicity above its injected value.
Corrected mass/energy errors remain below 2e-11. Changing surface tension rejects
restart before history output. Existing application and MPI API conservation,
rollback, migration and checkpoint cases also pass. Provenance, contract, Python
syntax and whitespace checks pass.

The application case demonstrates ETAB deformation-state continuation, not an
ETAB fracture in this very short interval; active fracture and component-mass
preservation were checked in the preceding MPI API increment. No selected-backend,
full matrix, physical injector validation or production scaling is claimed.
