# MPI AMR physical-boundary application and restart — 2026-09-17

Source: parent `6de364b77d49624d6cf6ca450ded88a7249be1d4` plus this change.
Logs, source hashes and the run artifact path are under
`build/amr-mpi-physical-app-checks/`. Reused GNU Fortran 13.3.0 Debug / Open MPI
4.1.6 configuration, with `/usr/bin/mpiexec.openmpi`.

The fixed and selected MPI AMR frontends now accept physical boundary inputs.
The shared application checks boundary settings collectively, constructs inlet
states, applies inlet acoustic/diffusion and wall thermal timestep limits, and
forwards physical data to the sparse full split. Accepted global boundary
transfers accumulate once and enter the conservation diagnostics. Existing
physical checkpoint schemas persist this ledger; root restores it and broadcasts
it before stepping. The checkpoint boundary signature is now a public shared
routine, without changing its serialized representation.

The wall thermal timestep helper now accepts one-column slabs. It visits each
face independently so both walls contribute when a cell touches both ends.
A focused unit fixture checks that two equal wall contributions halve the limit.

The application regression compares serial reference output with a four-rank
checkpoint followed by a two-rank restart. Wall, inflow and dynamic patch cases
all produce byte-identical coarse and fine CSVs. Changed wall/inlet settings
are rejected before output. Dynamic regridding changes the interior patch before
the checkpoint. The boundary ledger is compared with relative tolerance 2e-11
and absolute floor 1e-25; momentum components share their vector's maximum
component scale. An initial new-test failure used individual net-component
scales: cancelling pressure impulses produced a 2.9e-24 x-momentum difference
despite exact state output. That failed run is retained in
`ctest-momentum-scale-failure.log`; no existing test tolerance was changed.

Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target pelef_mpi_amr_reactive_3d pelef_amr_reactive_3d test_amr_reactive_transport_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^regression_mpi_amr_physical_boundary_restart_3d$' --output-on-failure
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^unit_amr_reactive_transport_3d$' --output-on-failure
python3 -m py_compile tools/check_amr_physical_boundary_restart.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Result: 2/2 focused tests passed, no skips (application 52.39 s, transport unit
17.08 s). Provenance, Python syntax, project contract and diff checks passed.
Run artifacts:
`build/spray-mpi-migration-gnu-debug/tests/mpi_amr_physical_boundary_restart_3d/run-l5n007u7/`.

This verifies the fixed-mechanism application with chemistry disabled and a
strictly interior two-level patch. Selected-mechanism application qualification,
boundary-touching fine patches, general hierarchies and physical/scaling
validation are outside this check.
