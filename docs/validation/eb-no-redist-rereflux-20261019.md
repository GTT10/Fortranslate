# Ordinary held EB MOL flux re-reflux, 2026-10-19

The nonuniform reacting EB-AMR discrepancy includes a missing operation in
PeleF's ordinary held-energy closure: pinned AMReX redistributes a cut cell's
coarse/fine flux mismatch even when PeleC selects `NoRedist` for the spatial
operator. This increment connects that operation for ordinary paired gas/UEINT
closures in serial and MPI. It does not claim complete trajectory parity.

## Causal diagnosis

Artifacts are under `build/eb-mol-spatial-diagnosis-20261019/`.
A diagnostic copy of pinned `Diffusion.cpp` exports the initial spatial
operator; no pinned source or reference executable is modified. The initial
coarse density, momentum and total energy match the imported Fortran state
exactly. The initial spatial-rate differences are much smaller than the
previous coarse-cut-cell endpoint discrepancy. The largest endpoint errors
are at coarse/fine interfaces, including their periodic images.

Pinned `AMReX_EBFluxRegister.cpp::Reflux` calls `eb_rereflux_from_crse`
independently of the spatial redistribution setting. For a cut boundary cell,
it retains kappa times the density increment and distributes the remaining
conserved inventory over connected neighbors, including fine-covered coarse
recipients. PeleF already implements this kernel and its distributed routing
for spatial ledgers; its ordinary flux-only closure had applied the raw
conservative mismatch directly.

A second diagnostic executable bypasses only this native re-reflux operation.
The transport-off, one-step comparison with the previous Fortran trajectory
then changes as follows:

| Maximum composite difference | Pinned native | Diagnostic without re-reflux |
| --- | ---: | ---: |
| Density, kg/m3 | 5.28760e-6 | 2.30101e-7 |
| x momentum, kg/m2/s | 2.46586e-3 | 5.07347e-5 |
| y momentum, kg/m2/s | 2.45546e-3 | 2.30162e-4 |
| Temperature, K | 1.48476e-3 | 1.48476e-3 |

The x-momentum discrepancy shrinks by about 49 times. Remaining fine-level
and temperature differences are independent work. The modified native binary
is a causal diagnostic, never the reference acceptance authority. Exact
compiler/link/run commands and its source copies are retained in the artifact
directory; `no-rereflux-comparison.json` records the numbers and worst cells.
Two initial diagnostic Fortran builds failed on missing module include paths;
the final build uses the original target's include arguments.

## Implementation boundary

`initialize_flux_rereflux_register_3d` builds a private interface ledger from
an already integrated density mismatch. It accepts only coarse face-boundary
cells; covered coarse cells and ordinary same-level cells contribute nothing.
Invalid source data leave the caller's ledger unchanged.

Ordinary held closure first assembles the conservative gas and UEINT mismatch,
then packs both into one ledger and applies the existing connectivity/periodic
re-reflux assembly. The difference between redistributed and raw increments
is delivered through the existing coarse/fine route, including fine-covered
recipients and ancestor spills. Both quantities undergo the same EOS admission
before publication. A failed admission retains the existing conservative
repair and transactional retry. MPI uses owner-local registers and collective
failure agreement, without gathering the gas tree.

This change is scoped to the ordinary paired held-energy closure. Spatial
StateRedist ledgers and unpaired/gas-only closure paths are not newly qualified
by this result; the complete flux-plus-spatial-register contract still needs
its own audit. Native input contracts, source SHA pins and acceptance limits
are unchanged.

## Verification

Working source is dirty `main` based on
`0a0e00ebe9138145fbb0e0e745665d0d768df009`. `before-identity.json`, `before/`
and `increment.patch` isolate this increment from earlier work. Toolchain is
GNU 13.3 Debug with runtime checks, OpenMPI 4.1.6 and the existing SUNDIALS 7.2
Fortran build. Native diagnostics reuse the pinned Release build's objects
with separately compiled replacement objects and separately named binaries.

Completed: LiDryer provenance/regeneration preflight, Python syntax checks,
project contract check and `git diff --check`. The serial reference probe
build completed. The final serial/MPI application and unit builds complete. Two focused CTests
pass: `unit_mpi_eb_interface_source_3d_3` (0.71 s) and
`unit_mpi_eb_native_mol_partial_recovery_3d_3` (12.82 s). The latter retains
source-frozen recovery, serial/MPI agreement, independent inventories and
provider-refusal rollback. No tests are skipped.

The first interface-unit attempt failed because its raw mismatch adapter
accepted diagonal cells from the redistribution mask. The adapter now also
requires a shared face with a refined cell; the targeted test covers this
exclusion and malformed-source rollback. Both reference and application checks
were rerun after that numerical correction. The preceding reference/application
runs are retained separately and are not additional final-pass counts.

The final hot-CVODE recursive application harness passes five successful runs
(serial, reaction off, MPI 2, stop on MPI 2, resume on MPI 3) and two expected
input refusals. Intermediate and final gas snapshots match byte-for-byte across
serial/MPI/changed-rank continuation. Independent mass/species/element/energy
audits pass. Its record is
`final-application-runs/reacting_restart_wd5tcb0p/checks.json`.

The final three-step reference probe completes but still fails 21 field gates.

The unchanged native capture is reused from
`build/tgreact-eb-amr-reference-20261019/runs/matched_eb_amr_62p3fjsx/reference`.
The new result is `final-reference-runs/matched_eb_amr_sjpmrwth/checks.json` under
the artifact directory. At 300 ns:

| Maximum field difference | Before | After |
| --- | ---: | ---: |
| Density, kg/m3 | 1.29714e-5 | 1.81190e-7 |
| x momentum, kg/m2/s | 6.62163e-3 | 6.55983e-4 |
| Temperature, K | 1.49761e-2 | 2.77874e-3 |

All three worst cells are now on the fine level. Species and H2-consumption
gates still pass, as do mass/element and numerical-energy budget checks.
`reference-summary.json` records the maxima across all four output times.
Further work should isolate fine-level halo interpolation and held-energy/EOS
conversion; conservation alone does not qualify the remaining field mismatch.

Exact source SHA256 entries are in `workspace-sha256.json`; touched numerical
sources and the unit fixture are also copied under `source/`. Final commands,
result paths and executable SHA256 values are recorded in `verification.json`
and each harness report. The native diagnostics are separate from the pinned
capture reused by the unchanged acceptance checker.
