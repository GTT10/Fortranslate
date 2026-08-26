# Decision 0135: route compact sparse child transport context

## Context

Owner-tiled sparse transport already limits root numerical work to finite row
bands and stores flux mismatch on patch support. After those stages, however,
each distinct fine-child owner still received complete root start/end states,
temperatures, and x/y fluxes. The child reads only four coarse exterior edges
for boundary reconstruction and the coarse/fine interface fluxes for its
register. It needed a complete corrected-root workspace only because reflux
ran on the child.

## Decision

Represent boundary input as raw start/end conserved state and temperature
samples on the four fine-patch edges. Preserve the established interpolation
and EOS recovery by reconstructing the exterior from those raw samples at each
substep. Initialize the compact flux register and accumulate coarse interface
fluxes on the root physics owner before sending the context.

The child owner advances the fine patch and accumulates fine interface fluxes,
then returns its fine state, temperature, and register correction. Apply
reactive reflux on the root owner in child order and return only the corrected
fine state and temperature. Count each packed payload as one transfer.

## Consequences

Remote transport child owners allocate no complete root state, temperature,
or flux field. The former bundle broadcast and patch-plus-two correction
round trips leave the transport path. A remote child now uses three messages
per Euler stage, and the context payload is bounded by child perimeter plus
patch-local register support. Root-owner reflux retains the established
deterministic correction order and complete-root compatibility kernel. Hydro
is unchanged.
