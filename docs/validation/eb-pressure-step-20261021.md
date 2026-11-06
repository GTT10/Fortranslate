# Paired pressure: one hydro increment and reflux, 2026-10-21

The diagnostic paired spatial-pressure reconstruction now advances a
forward-Euler hydro-only increment and applies the production gas and scalar
coarse/fine reflux routines. It updates the original conservative state,
not the auxiliary pressure-reconstructed energy arrays. Gas synchronization
also conservatively averages covered coarse cells; held energy is synchronized
through its independent scalar fluxes.

This is a one-step integration experiment following the
[initial spatial audit](eb-pressure-ghost-20261021.md). It is not the production
SDC integrator, recursive subcycling, a reacting trajectory, or MPI adoption.
No chemistry/transport source or energy reset is applied. The diagnostic
performs its updates on private copies and does not publish them to the
application state.

## Controls and accounting

The same planar 8-grid stationary import is used at dt=50 ns. Controls add
constant velocities (10,0,0), (-10,0,0) and (0,10,0) m/s, adding the matching
kinetic energy to total energy and keeping held thermal energy unchanged.
The x/y directions are periodic and the z boundaries are stationary slip
walls. Each control executes both the unchanged and paired spatial policy.
The moving controls must show a nonzero state increment.

Only composite leaves contribute to inventories: covered coarse cells are
excluded. The Fortran diagnostic records every conservative component and
held energy before/after reflux. The independent Python checker verifies
mass, total energy, each species, held energy and the species sum. Momentum
is not asserted to have a zero external budget because of wall forces.
The maximum E-minus-held-minus-kinetic difference and total/held EOS pressure
difference are recorded separately from integral conservation.

Acceptance limits are 1e-12 relative for integral inventories and species
closure, 1e-6 J/m3 for local E-H-KE consistency, and 1e-8 Pa for the stationary
paired pressure offset. Moving cases are not required to preserve pressure
exactly: conservative advection with the nonlinear EOS can create a finite
truncation error, which is reported against the original operator.

## Completed results

All four final controls pass the independent inventory checks. Maximum
conserved-component relative error is6.41e-15 and held-energy relative error
is2.29e-15 for the paired candidate. The moving increments are nonzero.

| Control | Original maximum pressure offset (Pa) | Paired offset (Pa) |
| --- | ---: | ---: |
| Stationary | 0.0196484 | 4.37e-11 |
| x velocity +10 m/s | 0.0196376 | 0.000327330 |
| x velocity -10 m/s | 0.0196346 | 0.000327330 |
| y velocity +10 m/s | 0.0197582 | 0.000104436 |

Paired local E-H-KE differences are at most5.83e-11 J/m3, and the two EOS
pressure reconstructions differ by at most4.37e-11 Pa. The original operator
also conserves the integral inventories; its different local E-H-KE behavior
is an observation, not a claim that all energy divergence must always vanish
in a general compressible flow.

Eight Python checks pass: four prior interface/restriction tests and four
independent budget checks (closure, energy leak, species/mass disagreement,
nonfinite rejection). Provenance, syntax, project contract and whitespace
checks pass. There are zero failures/skips in the final four-control run.

## Reproduction

The probe accepts `--audit-initial --separate-momentum --pressure-pair` and
writes `_baseline_pressure_step.csv`, `_paired_pressure_step.csv` and their
`_inventory.csv` companions. The existing spatial diagnostics remain available.

```sh
python3 tools/check_eb_pressure_step.py \
  --record build/eb-pressure-step-20261021/run-record.json \
  --output build/eb-pressure-step-20261021/checks.json
```

Commands, immutable copied executable, input identities, source hashes and
logs are retained in `build/eb-pressure-step-20261021/`. The earlier moving
pilot is separate from the final four-control record. GNU Fortran 13.3.0
Debug and the existing matching registration build are used; there was no
compiler/MPI configuration change.

## Remaining gates

A one-step conservative update does not prove multi-step stability or spatial
order. Next, connect the candidate through repeated staged hydro evaluation,
then check the actual held-MOL closure and pressure work, reaction/transport,
retry/rollback and serial/MPI equivalence. Keep the pinned reference default
and qualification distinct from this pressure-policy experiment.
