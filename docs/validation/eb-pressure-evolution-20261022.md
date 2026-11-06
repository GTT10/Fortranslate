# Paired pressure: staged hydro evolution, 2026-10-22

The pressure-policy probe now runs synchronous two-level SSPRK2. Both Euler
maps rebuild the production halos and donor pressures at their own stage
state, apply gas/scalar reflux, and then form the RK convex combination.
Temperature is recovered from held energy after stage updates and blending.
It does not reuse initial halos through an interval.

This extends the [one-step diagnostic](eb-pressure-step-20261021.md).
It is hydro-only: no reaction, transport, recursive subcycling, production
held-MOL reset or re-reflux/repair closure is invoked. The current experiment
uses the low-level conservative gas/scalar synchronization APIs. Its purpose
is to determine which remaining consistency condition blocks integration,
not to substitute SSPRK2 for the pinned SDC algorithm.

## Checks and scope

The stationary control and x=10 m/s advection run through200 ns. The moving
case also uses25 ns instead of50 ns at the same endpoint. This two-step-size
comparison measures sensitivity; it cannot establish temporal order.
At each accepted step, raw component and held-energy inventories are recorded
and independently checked. Missing endpoints, clock errors, duplicate or
missing inventory components, nonfinite outputs, species/mass disagreement
and integral leaks are rejected. Total energy is also decomposed into held,
kinetic and signed residual inventories; both signed and absolute local
residual integrals are retained to expose cancellation.

The earlier four-step moving pilot conserved inventories but accumulated a
local E-H-KE difference of1.55e-4 J/m3 at200 ns. Its maximum pressure offset
was0.001247 Pa. Therefore the one-step roundoff energy-gap result must not be
extrapolated into a multi-step energy-policy qualification.

The source of a required closure operation is explicit:
`close_held_prepared_eb_tree_3d` jointly routes/re-refluxes gas and held
increments, calls `reset_held_eb_tree_family_3d`, records the numerical
energy change, and then averages down. The current hydro-only probe bypasses
that operation. With the diagnostic's positive energies and eta2=2 policy,
the cell reset can replace total energy by held plus kinetic energy. Such a
change must be reported as an energy-policy correction, not silently counted
as conservative transport. No energy reset is added to hide the observed gap.

## Reproduction

The probe accepts `--audit-initial --separate-momentum --pressure-evolve`.
The V4 input step count and dt specify the duration. Both original and paired
policies start from the same imported state and use the same integrator.

```sh
python3 tools/check_eb_pressure_evolution.py \
  --prefix build/eb-pressure-evolve-20261021/moving50 \
  --input build/eb-pressure-evolve-20261021/moving50.dat \
  --output build/eb-pressure-evolve-20261021/moving50-checks.json
```

The work directory retains its20261021 start date across midnight:
`build/eb-pressure-evolve-20261021/`. The final frozen executable, commands,
input/binary hashes and reports belong to that directory. The earlier pilot
is separate from the final expanded energy-decomposition runs.

## Completed checks and pending comparison

The final stationary and moving50 ns runs both completed4 steps through200 ns.
Stationary paired pressure stays within4.37e-11 Pa; original pressure reaches
0.0613421 Pa. Moving paired pressure reaches0.00124674 Pa versus0.0612346 Pa
originally, with conserved-component relative error8.67e-15.

The moving paired local E-H-KE gap reaches1.55486e-4 J/m3. Its signed
composite integral is4.35e-20 J, while the absolute integral is2.13e-11 J.
Thus nearly cancelling local discrepancies are real despite the closed global
budget. The original absolute integral is4.26e-10 J. Do not report global
conservation as proof of local thermodynamic closure.

The25 ns comparison has completed and passed all inventory checks. At200 ns,
paired pressure changes from0.0012467403 to0.0012468065 Pa (+0.00531%);
maximum local E-H-KE gap changes from1.5548614e-4 to1.5461585e-4 J/m3
(-0.5597%); the absolute integral changes by-0.4547%. Thus halving dt does
not remove this discrepancy. These two sizes do not establish temporal order.
The next intervention is the existing cell dual-energy reset, with its signed
and absolute energy correction recorded separately from transport budgets.
Exact ratios are in `time-sensitivity.json`.
Five history-check tests and four inventory-check tests pass. Provenance,
Python syntax, project contract and whitespace checks pass.
