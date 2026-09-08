# Decision 0210: Make planar 3D EB restart transactional

## Context

The `0.216.0` executable can compose elementary chemistry, molecular
transport, stabilized planar EB hydro, and diagnostics, but every run starts
from the frozen density sheet.  A state-only dump would not be a reliable
restart: changing the NASA7 table, reaction mechanism, transport database,
geometry, or numerical controls could silently continue a state under a
different model.  Publishing state while a file is only partly parsed would
also leave the caller with an unusable half restart.

This increment is deliberately narrower than distributed or hierarchical
checkpoint I/O.  It owns one serial, single-level, exact-axis-plane EB state
and the execution history required by the existing public application.

## Decision

Add a self-identifying formatted checkpoint with an explicit schema and
terminal marker.  Store the complete NASA7 species records, elementary
reaction records, gas-transport records, immutable physics and numerical
configuration, planar EB geometry, conserved state, recovered temperature,
initial conserved/L1/element inventories, time and step, and cumulative
timestep/transport diagnostics.  Output paths, checkpoint scheduling, and
continuation limits are runtime policy rather than immutable physics; a
restart may extend `final_time` and `maximum_steps`, provided the stored clock
already fits the new limits.

Read every field into private candidates.  Validate the magic, schema,
dimensions, complete records, configuration and geometry fingerprints,
finite metadata, active-cell physical state, recovered NASA7 temperature,
and end marker before assigning any caller-owned field.  A missing,
truncated, incompatible, or physically inconsistent checkpoint therefore
leaves all state and history arguments bitwise unchanged.
Test each `iostat` result in a separate statement before inspecting the
corresponding record.  Fortran logical expressions are not required to
short-circuit, so combining an I/O failure and record validation could
otherwise inspect an undefined read target on early EOF.

Write a checkpoint only after a whole `R-T-H-T-R` step has committed.  The
application can stop immediately after that successful write or resume from
it and retain the pre-checkpoint minimum timestep, maximum diffusivity, and
minimum transport flux factor.  The public regression compares an
uninterrupted coupled run with a one-step checkpoint-stop and restart.

Normalize path components lexically before rejecting output/checkpoint and
output/restart collisions.  This covers repeated separators, `.` components,
and reducible `..` components without adding an operating-system dependency.

## Consequences

The supported serial planar EB workflow can now be split across processes
without changing its final deterministic CSV or cumulative diagnostic
history.  Incompatible thermodynamics, kinetics, transport, geometry, and
numerical settings fail before publication instead of being accepted on
species names alone.

The format is intentionally formatted and single-file.  It does not provide
parallel or scalable I/O, MPI ownership reconstruction, AMR topology,
regridding, general EB geometry, configurable external mechanisms, or forward
compatibility with a future schema.  Writes replace the configured file
directly; crash-atomic staging and rename are not claimed.  Those remain
separate workstreams.
Path validation is lexical: symlink, hardlink, bind-mount, case-insensitive
filesystem, and absolute-versus-current-directory aliases require filesystem
identity support and are not claimed by this portable input validator.
