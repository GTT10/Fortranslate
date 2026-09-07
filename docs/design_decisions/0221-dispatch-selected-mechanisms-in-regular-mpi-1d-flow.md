# 0221: Dispatch selected mechanisms in regular MPI 1D flow

## Status

Accepted for `0.228.0`.

## Context

The regular MPI 1D verification application already owns an uneven 19-cell
periodic decomposition, rank-local reactive state, communicator-wide stable
timestep and invariant reductions, transactional chemistry/transport/hydro
splitting, ordered gather, and deterministic CSV output. Its program unit and
MPI modules nevertheless imported the committed full-H2/O2 loaders and fixed
transport database type directly.

Linking a selected generated mechanism to `pelef_mpi_support` would also link
`pelef_core` and create duplicate module and symbol namespaces. Copying the
distributed time loop into a selected program would create a second MPI
algorithm. The selected bundle's declared chemistry policy must reach both
chemistry half-steps and fail collectively without publishing a rank-local
partial state.

## Decision

Extract `mpi_reactive_1d_application_mod`. Fixed and selected front ends load
and validate their own thermo, reaction, and transport arrays, then call this
one driver for deterministic initialization, uneven decomposition, adaptive
`R-T-H-T-R` advancement, collective invariants, gather, and dynamic-species
CSV output. Preserve the fixed command line as an optional output filename and
preserve its ten-species initialization and output byte-for-byte.

Make the MPI reactive and transport modules depend on
`gas_transport_species` rather than a fixed database loader. Thread an
optional `chemistry_integrator` through the distributed chemistry call,
Strang step, both chemistry half-steps, and adaptive retry. An absent argument
retains the fixed species-count policy. Require the resolved policy code to
agree across all ranks; an invalid or rank-disagreeing argument returns failure
on all ranks and leaves state and temperature unchanged.

Build `pelef_selected_mpi_reactive_1d_runtime` in a private Fortran module
directory on `pelef_selected_reactive_1d_runtime` and `MPI::MPI_Fortran`.
Configure and install `pelef_mpi_reactive_1d_selected` only when MPI is
enabled. Its link graph excludes `pelef_core` and committed mechanisms.

The selected verification case intentionally retains the fixed driver's
bundle-order composition wave and output-filename-only interface. Admit only
the exact ordered H2/H two-species profile or the pinned ordered ten-species
full-H2/O2 profile. The former supplies a nontrivial explicit-chemistry gate;
the latter retains the exact fixed initialization. This milestone does not
claim a general namelist-driven selected MPI application.

Require the MPI wrapper and launcher to resolve to one installation directory.
After linking the fixed and each selected coupled-MPI executable, inspect its
runtime dependency closure and reject unresolved libraries, conflicts, or
multiple versioned MPI ABI families.

## Consequences

Schema-1 bundles matching either qualified initialization profile can execute
the established regular MPI 1D verification path with generated explicit or
implicit chemistry and molecular transport. Fixture and pinned full-H2/O2
cases run with 1, 2, and 4 ranks. The full selected output must also equal the
fixed MPI output, while a separate nonzero adaptive-Strang unit requires exact
invalid/rank-disagreeing policy rollback.

This milestone does not add selected MPI AMR/EB, checkpoint/restart, general
input-driven initial or boundary conditions, runtime mechanism loading,
CVODE-backed CFD chemistry, mechanisms above 32 species, thread safety,
performance/scaling evidence, detailed fuels, external PeleC field parity, or
experiment validation.
