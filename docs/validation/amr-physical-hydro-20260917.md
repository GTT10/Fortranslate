# Serial two-level AMR physical hydro boundaries — 2026-09-17

Source: parent `42a737b` plus this change. Logs and source SHA256 hashes:
`build/amr-physical-hydro-checks/`. Reused GNU Fortran 13.3.0 Debug /
Open MPI 4.1.6 configuration; the exercised AMR API is serial.

`advance_amr_reactive_hydro_3d` accepts optional physical face definitions,
reservoir states and a conserved boundary-change vector (positive into the
hierarchy). The coarse boundary SSPRK2 solver optionally returns its mean
upper-face fluxes for reflux. Fine subcycling interpolates in time between
coarse states padded with physical ghost layers. Translating the interpolation
patch by two cells retains the existing interpolation kernel while preventing
periodic wrap at a physical coarse boundary. Reflux and average-down use the
original patch. Both levels and the boundary ledger commit only after success.
Omitting face metadata retains the existing periodic path.

The focused regression includes the existing periodic PCM/characteristic-PLM
checks and new physical-boundary cases for both reconstructions, HLLC and MC.
A ratio-two patch covers coarse indices 2:3 on a 12-cubed coarse grid. Slip-wall /
outflow and reservoir-inflow / outflow cases close composite conserved inventories
against the independently integrated boundary ledger within
`2e-11 * max(1, abs(initial_inventory))`. Doubling the far-outlet state leaves
all fine states and temperatures exactly unchanged. Negative inlet density
rejects the step and retains both input levels and temperatures exactly, with
zero returned boundary change.

Commands:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_amr_reactive_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^regression_amr_reactive_3d$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Result: focused regression passed 1/1 in 27.54 s; mechanism provenance, project
contract and `git diff --check` passed. No tests were skipped in this selection.

During development, a declaration-order compile error was corrected. The first
new rollback fixture incorrectly treated negative temperature guesses as invalid
conserved states; it now supplies negative density. The initial eight-cell
independence fixture allowed legitimate multistage PLM stencil propagation from
the far outlet. Enlarging that fixture separates stencil propagation from
periodic-wrap contamination without relaxing its exact-equality assertion.
Both failed test logs are retained. No Python source was changed.

Scope remains one strictly interior fine patch and serial hydro. Physical-boundary
chemistry/transport splitting, MPI and application input/restart integration,
fine patches touching domain faces, general multilevel/multipatch AMR and
physical validation are not established by this check.
