# Decision 0188: Run reacting full physics on boundary patch trees

## Context

The public arbitrary-depth x-upper cases qualify hydrodynamics and every
currently implemented dilute mixture-transport term through fresh and
restarted serial/sparse lifecycles. Chemistry remains disabled, so those cases
do not qualify reaction-modified states crossing the physical, recursive AMR,
process, and checkpoint boundaries together.

## Decision

Enable the elementary chemistry model in the existing fresh and split-run
boundary cases. Retain the full transport combination and require the complete
transactional `R-T-H-T-R` composition on all four populated levels. Reuse
serial and 1/2/4/8-rank fresh parity plus serial and cross-rank checkpoint
continuation rather than adding duplicate cases.

Do not change checkpoint schema 4. Its fingerprint already records chemistry
activation, model, relative and absolute tolerances, all transport controls,
and transport CFL.

## Consequences

The public x-upper lifecycle now exercises the selected elementary chemistry,
every currently implemented transport term, and hydro together while topology
is rebuilt, repartitioned, checkpointed, and restarted. The chemistry claim
does not imply arbitrary mechanism parsing or CVODE parity.
