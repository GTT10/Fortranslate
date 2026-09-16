# MPI physical-boundary AMR regridding — 2026-09-17

Source: parent `2b460ed` plus this change. Logs and source hashes:
`build/amr-mpi-boundary-regrid-checks/`.
Reused GNU Fortran 13.3.0 Debug / Open MPI 4.1.6 build with bounds checks.

MPI applications now accept `regrid_physical_boundaries = .true.`. The tagging
policy is checked collectively alongside the existing interval, buffer and
threshold. The root tagging/remap calls receive physical face and reservoir
data through the shared serial implementation. New patch geometry is broadcast,
sparse ownership is rebuilt, and remapped states are redistributed. The default
false flag retains the previous interior tagging/remap path. Existing dynamic
checkpoint signatures already retain the physical-tagging policy, so no schema
or serialized record layout changed.

The new regression compares uninterrupted serial execution against four-rank
checkpoint creation and two-rank restart after step-zero boundary expansion.
It includes a reacting full-H2/O2 isothermal moving-wall/transport case and an
elementary reacting corner case. Coarse/fine CSV bytes must match exactly;
MPI boundary ledgers retain the established tolerance. Removing the physical
tagging marker from a checkpoint must reject restart before output. Existing
boundary/chemistry mismatch and inactive-chemistry checks remain. The prior
fixed-chemistry restart regression also exercises the default interior dynamic
path in the same run.

Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target pelef_mpi_amr_reactive_3d pelef_amr_reactive_3d --parallel 4
python3 -m py_compile tools/check_amr_physical_boundary_restart.py
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^regression_mpi_amr_(boundary_touching_dynamic|fixed_chemistry)_restart_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Result: 2/2 tests passed, no failures or skips (existing fixed-chemistry
restart 65.77 s, new boundary-touching dynamic restart 55.53 s; total 121.33 s).
Provenance, Python syntax, project contract and diff checks passed. Final new
application artifacts:
`build/spray-mpi-migration-gnu-debug/tests/mpi_amr_boundary_touching_dynamic_restart_3d/run-5za8hvgs/`.

Regridding still gathers full levels on rank zero and is not a scalable
fully distributed remap. General multipatch/multilevel hierarchies and EB
coupling remain open. These checks are application/restart numerical evidence,
not external solver parity, scaling or physical validation. The selected
frontend shares the application path but was not independently qualified here.
