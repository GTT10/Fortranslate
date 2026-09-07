# Decision 0200: isolate static 3D AMR hydrodynamic synchronization

## Context

Milestone `0.206.0` qualified reacting hydro, chemistry, and molecular
transport on a serial periodic uniform 3D grid. Adding AMR, source splitting,
diffusive reflux, dynamic topology, restart, and MPI ownership together would
make a conservation failure difficult to attribute. The first 3D AMR boundary
therefore needs to isolate hyperbolic coarse/fine synchronization.

## Decision

Introduce one rectangular fine patch strictly inside the periodic coarse grid
with an integer refinement ratio. Initialize it independently or through PCM,
advance a provisional coarse SSPRK2 step, and retain the arithmetic mean of
the two stage fluxes in x, y, and z. Take `r` fine SSPRK2 substeps. Fill all
six fine boundaries from linearly time-interpolated provisional coarse states
and recover their NASA7 temperatures before each Riemann solve.

Average fine interface fluxes across `r^2` child faces and `r` substeps. Apply
their difference from the coarse time-averaged flux to the surrounding
uncovered coarse cells, then average the fine state onto covered parents.
Keep both levels in private candidates until reflux, restriction, and all EOS
recoveries succeed. Select a root interval no larger than the coarse stable
step or `r` times the fine stable step.

Qualify the implementation with geometry/unit gates, a nonuniform 3D
all-component conservation regression, uniform invariance, invalid-solver
rollback, and a public full-H2O2 coarse/fine CSV contract.

## Consequences

PeleF now has a runnable serial static two-level 3D hydrodynamic AMR boundary
with ratio subcycling, time-interpolated coarse ghost states, six-face reflux,
and conservative average-down. Chemistry and molecular transport are not yet
composed across levels. PCM spatial reconstruction/prolongation, a strictly
interior single patch, periodic coarse boundaries, serial ownership, and no
restart or EB remain explicit limitations.
