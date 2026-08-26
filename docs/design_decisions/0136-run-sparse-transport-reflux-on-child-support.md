# Decision 0136: run sparse transport reflux on child support

## Context

The compact child transport context removed complete root fields from remote
fine owners, but reactive reflux still ran on the root physics owner. Each
remote child therefore returned its complete evolved fine field with the flux
register and received the complete corrected fine field back. Reflux can alter
coarse state only on the patch-plus-one mismatch ring and one further
cardinal/diagonal recipient layer.

## Decision

Generalize nonreactive and reactive reflux to accept coarse arrays with global
lower bounds and any support containing the patch expanded by two coarse
cells. Preserve the complete-root entrypoints as wrappers. Restrict reactive
coarse temperature recovery to the supplied support while retaining complete
fine-patch recovery and transactional register reset.

Include the latest corrected patch-plus-two state and temperature with the
exterior/register context. Execute fine transport and reactive reflux on the
child owner. Keep corrected fine state there and return only corrected coarse
support to the root owner, which merges supports in child order.

## Consequences

A remote transport child uses two messages per Euler stage instead of three,
and fine state no longer crosses the root boundary for reflux. Payload storage
scales with child perimeter, compact register area, and patch-plus-two coarse
area. The root owner retains deterministic cumulative correction order and the
final full-root scatter. Hydro is unchanged.
