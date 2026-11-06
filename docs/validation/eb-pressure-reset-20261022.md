# Paired pressure with explicit dual-energy reset, 2026-10-22

The hydro-only SSPRK2 probe now has an explicit
`--pressure-evolve-reset` diagnostic. After each completed RK step it calls
the existing `reset_nasa7_dual_energy_cell` with the same positive-energy,
eta2=2 policy as the held reference probe. Fluid leaf states are reset on
private copies, then conservatively averaged down. Covered coarse states
are excluded from the numerical-energy inventory.

This tests the previously missing cell reset, not the complete prepared
held-MOL closure. The production re-reflux repair, reaction/transport,
recursive clocks, MPI and application checkpoint paths remain unmodified.
Both original and paired pressure policies receive the same reset schedule.
The earlier no-reset evidence is preserved in the
[evolution record](eb-pressure-evolution-20261022.md).

Every accepted history row records the cumulative signed energy correction,
cumulative sum of absolute cell corrections, raw total-energy change and
correction-accounted conservation residual. The independent checker subtracts
only the declared correction when testing total-energy balance and separately
checks raw total=held+kinetic+signed-gap decomposition. A reset is therefore
not silently presented as conservative transport. Held energy must remain
unchanged by this diagnostic reset; other policy branches are refused.

The execution chain runs the fixed binary and immediately performs the
independent history/budget check. It records command return codes and hashes.
The active tool call emits completion only after both stages finish, allowing
analysis to continue in the same assistant turn. This is not a configuration
change to Codex's user-facing turn-completion notification hook.

Artifacts: `build/eb-pressure-reset-20261022/` (frozen binary, driver,
record.json, CSVs, checker reports and logs).

## Completed results

Both stationary and moving inputs completed four50 ns steps through200 ns,
with both original and paired policies. All four histories pass the independent
correction-accounted gas/species/held inventory checks. Seven history tests
and four inventory tests pass; provenance, syntax, project contract and
whitespace checks pass. GNU Fortran13.3.0 Debug; no full configuration matrix,
MPI or native-reference rerun is claimed.

For the moving paired case, the accepted local E-H-KE gap becomes exactly
zero after every reset. The cumulative signed correction at200 ns is
7.42148e-18 J and the sum of absolute cell corrections is2.13014e-11 J.
The correction-accounted worst inventory error is8.67e-15. The original-policy
absolute correction is4.25550e-10 J, about20 times larger. These are explicitly
numerical energy corrections, not conservative transport fluxes.

The paired moving maximum pressure offset remains0.0012467403 Pa, unchanged
from the no-reset diagnostic at reported precision. Stationary paired pressure
remains at4.37e-11 Pa. Thus the accepted energy reset closes the local energy
identity but does not eliminate the residual advected-pressure truncation
error. It would be incorrect to attribute that pressure error solely to a
missing reset. Full prepared-MOL/re-reflux integration and spatial accuracy
remain open.

The tool completion event was received during the active turn and followed
by report analysis; the completed jobs were not handed back merely awaiting
manual confirmation. This does not establish post-turn automatic wakeup.
