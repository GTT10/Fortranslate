# Decision 0194: Freeze the completion contract and harden toolchains

## Context

The numerical suite had reached a broad serial and sparse-MPI 2D EB milestone,
but project-level completion was not mechanically defined. The CMake version
and runtime banner had drifted, the PeleC comparison branch was not frozen at
repository level, a mixed Intel/OpenMPI/GNU environment could pass MPI package
detection and fail later while reading `mpi_f08.mod`, and internal callback
trampolines caused several ELF targets to request an executable stack.

## Decision

Record the complete recursive PeleC reference snapshot in a machine-readable
manifest and make the completion workstreams and exit gates explicit. Add a
project-contract test that keeps CMake, runtime, README, validation, presets,
and the reference record aligned.

Require a usable `mpi_f08` module immediately after MPI discovery and report
the selected wrapper and compiler when they are incompatible. Pass explicit
geometry context through serial and sparse patch-tree regridding callbacks.
Move every callback whose address escapes a host procedure to module scope so
GNU Fortran does not emit stack trampolines. Inspect the final ELF program
headers in CTest and reject an executable `PT_GNU_STACK` segment.

## Consequences

MPI toolchain contamination now fails during configuration with an actionable
message. Geometry construction no longer depends on hidden host association,
which also makes callback state explicit and reentrant. Shipped serial/MPI
executables and the direct callback regression no longer request executable
stack permission. The new gates harden the existing subset; they do not add
3D, LES, detailed chemistry, particles, spray, or GPU capability.
