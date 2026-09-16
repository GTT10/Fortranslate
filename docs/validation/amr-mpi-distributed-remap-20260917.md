# Distributed sparse MPI AMR remap — 2026-09-17

Source: parent `57ccf01` plus this change. Logs and source hashes:
`build/amr-mpi-distributed-remap-checks/`.
Reused GNU Fortran 13.3.0 Debug / Open MPI 4.1.6 build with bounds checks.

The sparse remap API operates entirely on local coarse/fine arrays. It validates
old/new geometry and collective species/reconstruction metadata, restricts old
fine cells locally, exchanges two-layer coarse halos, and constructs children
for the new patch. Conserved-variable limited slopes, species density closure,
positivity backtracking and temperature recovery follow the serial remap.
Fine overlap is reused exactly. Its parent coarse cell remains on the same
rank because coarse ownership is fixed; unsupported coarse repartitioning is
rejected. Local final restriction and global composite-inventory reductions
retain the 5e-12 conservation guard. Candidate hierarchy, geometry and ownership
commit together only after all ranks succeed.

The application uses this API directly after distributed tagging. It no longer
gathers full fields or scatters remapped state at geometry changes. Fine
coordinate arrays are rebuilt locally. Initial/restart distribution and
checkpoint/CSV I/O retain their root aggregation; this is not a distributed-I/O
or load-balancing implementation.

A small four/eight-rank remap test compares all coarse/fine states and
temperatures exactly with serial remapping. It exercises movement, growth,
shrinkage, lower/upper physical corners and full-domain coverage at refinement
ratios two and three, PCM/PLM, fine-scale perturbations and ranks gaining or
losing fine ownership. Invalid or rank-mismatched target geometry and a
one-rank invalid inlet preserve both levels, geometry and ownership. Existing
fixed-chemistry and boundary-touching dynamic application regressions compare
serial CSVs with four-to-two-rank checkpoint/restart, including their context
rejection probes.

Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target test_mpi_amr_regrid_3d pelef_mpi_amr_reactive_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^(mpi_amr_regrid_3d_np(4|8)|regression_mpi_amr_(boundary_touching_dynamic|fixed_chemistry)_restart_3d)$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Result: 4/4 tests passed, no failures or skips (remap four ranks 1.17 s,
eight ranks 1.07 s; fixed-chemistry restart 67.89 s; boundary dynamic restart
56.86 s; total 127.01 s). Provenance, project contract and diff checks passed.
No Python source changed. The build succeeded with real-equality warnings.
Final application artifacts:
`build/spray-mpi-migration-gnu-debug/tests/mpi_amr_fixed_chemistry_restart_3d/run-nglbb7dt/`.
`build/spray-mpi-migration-gnu-debug/tests/mpi_amr_boundary_touching_dynamic_restart_3d/run-_9hklyt0/`.

General patches/levels, coarse repartitioning, distributed I/O and EB coupling
remain open. These checks establish numerical equivalence and transactional
behavior, not measured scaling, external solver parity or physical validation.
No full configuration matrix was run.
