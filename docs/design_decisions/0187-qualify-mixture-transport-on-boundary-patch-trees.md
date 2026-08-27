# Decision 0187: Qualify mixture transport on boundary patch trees

## Context

The public arbitrary-depth x-upper cases now qualify recursive Fourier
conduction through fresh and restarted serial/sparse lifecycles. The same
transport kernel already supports viscosity, mixture-averaged species
diffusion, barodiffusion, correction velocity, and species enthalpy flux on
interior trees, but those terms are not exercised at the physical-side path.

## Decision

Enable every currently implemented molecular-transport term in the existing
fresh and split-run boundary cases. Reuse the established recursive SSPRK2
schedule, current-fine physical exterior, coarse-time context elsewhere,
`r^2` child subcycling, diffusive registers, reflux, and average-down. Require
serial and 1/2/4/8-rank fresh parity plus serial and cross-rank checkpoint
continuation.

Do not add duplicate cases or change checkpoint schema 4. Its fingerprint
already records viscosity, thermal conduction, species diffusion,
barodiffusion, and transport CFL.

## Consequences

The public boundary lifecycle now exercises the complete currently
implemented dilute mixture-transport combination across physical, recursive
AMR, process, and restart boundaries. This is not a claim of Stefan--Maxwell,
Soret, Dufour, or complete PelePhysics transport parity.
