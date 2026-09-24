# Regular MPI spray hydro PPM operator — 2026-09-25

Source: base `67e0850`, implementation `fcd03e2`. GNU Fortran 13.3.0,
Open MPI 4.1.6 (`/usr/bin/mpiexec`), Ninja Debug, Python 3.12.3. Reused fixed
MPI tests-enabled `build/spray-amr-reactive-debug`. Artifacts:
`build/ppm-regular-mpi-hydro-checks/`.

`advance_mpi_spray_hydro` accepts an optional trailing `reconstruction`.
The default remains `characteristic_plm`; `characteristic_ppm` selects one CTU
step. All ranks agree on the selector before exchange. PPM uses the existing
variable-width MPI LES boundary sampler to obtain four inviscid ghost layers,
including when the slab is narrower than four cells. It computes owned-cell
fluxes and EOS recovery privately and reduces physical surface fluxes for the
boundary ledger. Failures retain U/T and leave the optional ledger zero.
Absent face arguments retain periodic behavior. Particle coupling and the
regular MPI application do not yet forward/expose this selector.

## Completed checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_ppm --parallel 2`: passed.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset,
  `HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`:
  `ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_ppm_4$' --output-on-failure`:
  **1/1 passed**, 2.12 s test / 2.16 s total, four ranks.
- An 8x4x4 nonuniform-density H2/O2/N2 field yields two x cells per rank.
  Periodic operation without explicit faces, outflow/slip-wall faces, low-x
  Mach inlet and high-z Mach inlet are tested. Owned U/T agree exactly with
  the regular serial physical-boundary PPM operator. Boundary ledgers agree
  within 2e-13 scaled tolerance. Oversized dt and one-rank PLM/PPM disagreement
  reject collectively with unchanged U/T; the failed dt ledger remains zero.
- `python3 -m py_compile tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, `git diff --check`: passed.

No failed tests occurred; compiler warnings remain. Selected builds, particle
coupling, application/restart, molecular transport/LES/chemistry in this fixture,
full matrix, convergence and physical validation were not run here. This does
not establish exact reference primitive-pressure parity or spray validation.
