# Serial AMR physical-boundary molecular transport — 2026-09-17

Source: parent `9e3a0b3` plus this change. Source hashes and logs:
`build/amr-physical-transport-checks/`. Existing GNU Fortran 13.3.0 Debug /
Open MPI 4.1.6 configuration; serial fixed H2/O2 API tested here.

The two-level transport API accepts optional face definitions, reservoir states,
wall temperature/velocity and a conserved boundary-change vector. The coarse
Euler stage uses existing physical molecular fluxes, retaining both domain-edge
faces for the update and surface budget and upper-cell faces for reflux.
Fine transport subcycles interpolate between physically padded coarse states;
the translated patch prevents periodic wrapping in near-wall interpolation.
Existing fine-interface limiting, reflux and average-down retain their original
geometry. SSPRK2 returns the mean of the two stage boundary transfers and
commits only after both stages and temperature recovery succeed.

The complete reaction / transport / hydro / transport / reaction API forwards
these settings and accumulates the two transport budgets plus the hydro budget.
It publishes states and budget only after the entire split succeeds. Omitting
face metadata retains periodic stepping; disabled transport forwards directly
to the reaction/hydro split with the same boundary options.

The focused unit test retains its periodic, per-process, six-face adversarial
interface-limiter and ratio-three cases. New cases on a 4-cubed coarse grid with
a ratio-two interior patch exercise a 600 K isothermal wall moving tangentially
at 20 m/s, opposite outflow, and a separate prescribed reservoir inlet. Viscosity,
conduction, species diffusion and barodiffusion are enabled. Checks cover nonzero
boundary energy transfer, cooling direction, composite conserved inventory
against boundary flux, limiter bounds and exact rollback with zero diagnostics
on invalid wall temperature or inlet density. Full split cases exercise both
chemistry disabled and enabled, including H/O/N boundary balance and rollback
of the first chemistry half when transport rejects the wall input.

Conserved balance uses relative 2e-11 against initial inventory or boundary
transfer (floor 1e-20); momentum uses a shared vector-component scale. Reacting
species are assessed through element balance, relative 2e-10 against initial
elements (floor 1e-30). No existing tolerance was changed.

Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target test_amr_reactive_transport_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^unit_amr_reactive_transport_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Result: the focused unit test passed 1/1 in 17.37 s, with no skips. Mechanism
provenance, project contract and `git diff --check` passed.

The first new inlet fixture incorrectly retained a nonzero wall velocity on an
inflow face. The input validator rejected it; the corrected fixture sets wall
velocity only on the wall. `diagnostic.log` retains that failure. Temporary
stage diagnostics were removed. No Python source changed.

This establishes bounded serial numerical behavior, not experimental validation.
Physical-boundary MPI stepping and application inputs/restart, patches touching
domain faces and general multilevel/multipatch AMR remain outside this scope.
