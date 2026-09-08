# 0213: Run selected normalized mechanisms in a generic 0D reactor

## Status

Accepted for `0.220.0`.

## Context

Milestone `0.219.0` proved that a normalized mechanism selected at CMake
configuration time could be source-checked, generated, compiled, installed,
and exercised through a probe. It did not advance a physical state. The fixed
full-H2/O2 application already had an adaptive implicit constant-volume
integrator, but directly named committed thermo and mechanism modules.

Linking an arbitrary selected module to `pelef_core` is unsafe when the bundle
reuses a committed generated module name. The application also needs a stable
configuration contract independent of the bundle's species order.

## Decision

When `PELEF_MECHANISM_BUNDLE` is set, configure and install
`pelef0d_selected` alongside `pelef_mechanism_probe`. Link the generated
mechanism to a small `pelef_selected_runtime_core` containing only generic
precision, constants, NASA7 mixture properties, elementary kinetics, the
constant-volume solver, and selected-reactor configuration. Give both runtime
and generated targets distinct Fortran module directories; do not link the
committed generated mechanisms into this executable.

Read initial composition as exact species names paired with mole fractions,
then map it into generated bundle order. Reject nonfinite integration values,
invalid or nonzero unused array tails, duplicate or unknown names, invalid
reaction records, molecular-mass imbalance, and lack of a common NASA7
temperature interval before opening output. Reject portable lexical
input/output aliases after reducing `.`, repeated separators, and reducible
`..` components. Validate reaction element and molecular-mass balance from
normalized thermo records during generation as a second boundary.

Require `initial_time_step` to lie within the configured minimum and maximum
and require `output_interval` to be at least the minimum. Interpret
`minimum_time_step` as the lower bound for an adaptive reduction. A final or
output scheduler fragment may be shorter than that bound only when the solver
accepts the complete requested fragment without reducing it further.

Keep the generator's legacy kinetics-only input mode, but do not admit those
inputs to the CMake build helper. The helper and both configured executables
require schema 1, explicit interface names, and one NASA7 thermo and primitive
transport record per species, with an explicit configure-time diagnostic for
the legacy boundary.

Reuse `advance_constant_volume_implicit_adaptive` with the dynamically loaded
species and reaction arrays. Use the generated production-rate kernel only for
the dynamic CSV molar-rate columns, whose `wdot_*` values are in
`kmol/(m^3 s)`; do not describe the integrator as a generated callback path.

Require a reactive two-species end-to-end case with independent conservation,
activity, exact name-mapping, startup-failure, and repeat-run byte-identity
checks. Require the selected pinned full-H2/O2 executable to reproduce the
fixed full-H2/O2 CSV byte-for-byte. Install and non-executable-stack audit both
selected executables.

## Consequences

A valid configure-time selected schema-1 mechanism with 2--32 species can now
run as a serial adiabatic constant-volume reactor without editing fixed source
modules. A selected full-H2/O2 build may safely declare the same module names
as the committed adapters because no executable links both definitions.

The application is not runtime YAML/JSON loading, a mechanism switch inside
one binary, CFD/transport/AMR/EB/MPI dispatch, CVODE/SUNDIALS parity,
production-performance qualification, or physical validation of arbitrary or
detailed-fuel chemistry. Output publication after a later integration failure
has the same direct-CSV boundary as the existing fixed reactor applications.
Path identity is lexical; symlink, hardlink, mount, case-folding, and
absolute-versus-current-directory aliases are not claimed.
