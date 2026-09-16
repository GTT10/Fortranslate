# MPI Smagorinsky coupling and restart — 2026-09-17

Source: parent `fb979bb` plus this change. Source hashes and build/preflight/test
logs in `build/spray-mpi-les-checks/`. Existing GNU Fortran 13.3.0 Debug / Open MPI
4.1.6 configuration.

A bounded two-halo wrapper now exposes the existing serial SGS face kernel.
The MPI periodic x-slab updater exchanges gas/temperature halos, computes the
same deviatoric stress, heat and species fluxes, and advances SSPRK2 with collective
thermodynamic acceptance. Each stage rejects diffusion number above 0.2 before
committing state. The coupled driver inserts SGS half-steps around hydro, inside
molecular-transport half-steps. The application's existing whole-step retry handles
rejections, including transactional particle/ID state.

`les_cs`, `les_prt`, `les_sct` configure the MPI application. Default-disabled
contexts are preserved, and nondefault coefficients enter checkpoint identity.
No new persistent SGS field is introduced by this algebraic model.

Commands:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi pelef_mpi_spray_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_(application|migration)_[24]$' --output-on-failure
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_application_[24]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
python3 -m py_compile tests/spray/test_mpi_spray_application.py
git diff --check
```

API tests compare nonuniform velocity/temperature/composition fields with the
serial SGS SSPRK2 reference on uneven 2/4-rank partitions (2e-12 relative field
tolerance, 1e-8 K temperature). The reference conserves total gas inventories;
MPI results agree with it and visibly change the gas. An excessive diffusion
step rejects without changing gas or temperature. Both MPI API entries passed.

The first application check failed its added temperature-effect assertion:
10 m/s injection over 2e-8 s gave exactly zero resolved temperature difference
and a maximum absolute x-momentum difference 5.90e-17 from the SGS-disabled run.
The implementation and changed-rank comparison passed up to that assertion.
The retained log is `ctest-initial.log`; the initial four-entry run had two API
passes and two application failures (72.27 s). The effect test was changed to
1000 m/s injection with its own identical SGS-disabled baseline and an explicit
x-momentum effect check. Conservation and restart tolerances were unchanged.

The application test also covers 2↔4-rank continuation from a 1e-8 s checkpoint
through 2e-8 s (relative 2e-11 / absolute 1e-25), mass/energy balance below 2e-11,
and rejection of changed Cs before output history. These are numerical coupling
checks, not a validated LES, wall model, turbulent dispersion, AMR/EB or production
scaling qualification. No full matrix or selected-mechanism SGS rerun was made.

The revised 2-rank effect case gives maximum absolute x-momentum change
6.80e-11. An initial compile-time control-array length mismatch was corrected
before numerical tests; no model coefficients or numerical tolerances were changed.

Final evidence: MPI API entries 2/2 passed (4.76/4.95 s); revised application
entries 2/2 passed (65.43 s total), with no final failures/skips. Provenance,
project contract, Python syntax and whitespace checks passed.
