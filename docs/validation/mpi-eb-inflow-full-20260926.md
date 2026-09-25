# MPI EB prescribed inflow transport and full split — 2026-09-26

## Behavior

Source 20b1373, parent 2c8685f, local main.

Common and recursive MPI EB molecular transport and full reacting-split APIs
now forward explicit inflow/outflow labels and prescribed conserved states.
Explicit open faces participate in transport and exterior-flux accounting;
legacy absent-label behavior remains unchanged. MPI timestep selection also
forwards these controls into the existing donor speed/diffusivity bounds.
Any rank supplying explicit data triggers collective halo validation.

Both transport halves and the hydro update use the same input controls.
The existing candidate state transaction preserves gas, temperature, limiter
and all ledgers if a later stage fails. Inlet temperatures remain EOS guesses;
composition comes from the prescribed conserved states.

## Verification

GNU Fortran 13.3.0, Open MPI 4.1.6 via /usr/bin/mpiexec, Ninja Debug,
Python 3.12.3, existing tests-enabled build/spray-amr-reactive-debug.
No compiler/backend changes.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_subcycle_full_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_ctu-inflow(-common)?_full_3d_3$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Each fixture uses five patches on three levels, CTU with nonzero artificial
viscosity, molecular transport and reacting gas. Explicit chemistry uses
distributed ownership; implicit chemistry puts all patches on the last rank,
leaving two idle ranks. State, temperature, boundary ledger and limiter must
match serial execution exactly. Wall/chemical ledger and total-inventory
residuals use the existing 3e-14 tolerance. MPI/serial inlet timestep bounds
are compared at 1e-14 relative tolerance.

Failure checks include invalid hydro after private reaction/transport stages,
one-rank nonfinite prescribed inlet temperature after a private reaction half,
and mismatched coupling controls. State, temperature and ledgers must remain
unchanged on failure. No numerical tolerances or golden data are changed.

Artifacts: build/mpi-eb-inflow-full-checks/build.log and LastTest.log.
Mechanism SHA: 0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4.

## Results

**2/2 passed, 0 failed, 0 skipped**, 196.47 s total: recursive 109.83 s,
common 86.59 s. Both explicit and implicit chemistry fixtures passed serial
parity, timestep and conservation checks. Late-stage and collective-input
failure cases preserved state and ledgers. No failed build or test attempts.
Provenance, Python syntax, project contract and whitespace checks passed.
Existing compiler warnings remain.

## Limits

These are fixed-mechanism MPI library checks. Tree application input,
checkpoint identity and spray wrapper forwarding remain to connect. No
periodic/wall tree dispatch, full matrix, installed audit, external solver
comparison or physical validation is established here. Full PeleC parity
remains unproven.
