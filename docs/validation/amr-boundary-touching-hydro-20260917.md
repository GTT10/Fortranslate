# Boundary-touching two-level AMR hydro — 2026-09-17

Source: parent `48ab976` plus this change. Logs and source hashes:
`build/amr-boundary-hydro-checks/`. Reused GNU Fortran 13.3.0 Debug build
with bounds and floating-point checks.

The serial two-level hydro API now accepts patches touching physical domain
faces. Each fine SSPRK2 stage maps physical ghost coordinates into the fine
state, applying reservoir inflow, outflow or normal wall reflection. Coordinates
outside an interior coarse/fine interface still use time-interpolated coarse
PCM/PLM data. This includes corners shared by physical and coarse/fine faces.
Hydro keeps the established normal-reflection rule for slip/no-slip/isothermal
walls; tangential viscous wall conditions belong to the transport operator.

The coarse physical SSPRK2 routine can now return all face fluxes including
lower domain faces. Its existing upper-face output shape remains supported.
Fine subcycle fluxes feed the boundary-aware reflux operator, replacing the
covered coarse surface contribution in the composite boundary ledger. Average
down and temperature recovery finish the candidate before either level or its
diagnostics commits. Reflux magnitude is now also private until final success.
Periodic seams still require separate ownership handling and are rejected when
touched by the patch. Interior patch arithmetic follows the existing path.

The focused hierarchy fixture covers all six faces, a full-domain patch and a
three-face corner; it uses PCM and characteristic PLM, ratios two and three,
nonuniform fine density, nonzero three-component velocity, inflow/outflow and
reflecting walls. It checks conserved composite inventories against the returned
surface transfer. Full-domain refinement is compared to repeated regular-grid
physical SSPRK2 steps at the fine resolution, including exact fine states and
temperatures. Invalid inlet data and unsupported periodic seams must preserve
both levels and zero diagnostics. Existing interior chemistry/split tests run
in the same executable.

The first new fixture run failed after its rejection probe restored inlet
density from a covered coarse cell rather than the saved reservoir value.
The original inlet density is now retained explicitly; the failed log is kept
as `ctest-inlet-fixture-failure.log`. Solver tolerances were not changed.

Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target test_amr_reactive_chemistry_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^unit_amr_reactive_chemistry_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Result: 1/1 test passed with no skips (36.25 s; CTest total 36.28 s).
Provenance, project contract and diff checks passed. No Python source changed.

Boundary-touching reaction/transport splits, MPI, regridding, checkpoints and
application inputs remain to be connected. This is numerical hydro evidence,
not general AMR or physical validation.
