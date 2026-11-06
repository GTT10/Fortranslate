# Held EB MOL coarse-to-fine halo interpolation, 2026-10-19

The remaining fine-level discrepancy after the ordinary flux re-reflux fix
includes a missing coarse-to-fine spatial interpolation. The held-energy MOL
path now uses a component-shared linear limiter where the coarse cell has all
27 EB connections. This is implemented for serial and owner-local MPI gas
storage. It is not a claim of complete EB-AMR trajectory parity.

## Diagnosis and numerical boundary

Pinned PeleC selects `eb_mf_lincc_interp` with its default state interpolation
settings. `AMReX_EBMFInterpolater.cpp` first applies the common-component
limited linear interpolater and restores PCM when the coarse cell has fewer
than 27 connection bits. PeleF previously copied the coarse cell everywhere,
although it already interpolated the ancestor endpoints in time.

`plan_eb_linear_halo_3d` uses the donor's extended geometry and the existing
open-face connectivity algorithm. It resolves axial state donors across
same-level peers and periodic seams, including reflected momentum at walls.
The actual gas, temperature and stored internal energy share one limiter per
direction, after ancestor time interpolation. A blocked EB stencil stays PCM.
Inlet-adjacent stencils and missing same-level donor neighborhoods also retain
the previous PCM behavior: recursive spatial stencils and `ext_dir` slopes
remain future work. Gas-only and non-MOL halo callers are unchanged.

Serial and MPI use the same planning and interpolation routines. MPI fetches
only the six requested axial cells, preserves owner-local gas storage and
collectively rejects failed planning or interpolation before publishing the
candidate halos. Held energy remains an independently transported quantity;
this change does not replace it with total energy minus kinetic energy.

A copied native `Diffusion.cpp` exports initial coarse/fine input halos before
the spatial operator. The pinned source, normal executable and acceptance
limits are unchanged. Fine valid gas and held energy agree exactly before
this change, while coarse-sourced fine halos differ. With the new spatial
interpolation, the first fine spatial evaluation gives:

| Maximum absolute initial fine-halo difference | Previous PCM | Linear |
| --- | ---: | ---: |
| Density, kg/m3 | 1.278633e-4 | 0 |
| x momentum, kg/m2/s | 1.206166e-2 | 1.110223e-16 |
| y momentum, kg/m2/s | 2.710282e-2 | 1.665335e-16 |
| Total energy, J/m3 | 6.990449e1 | 2.910383e-11 |
| Stored internal energy, J/m3 | 6.985479e1 | 0 |
| Temperature, K | 1.636051 | 1.689758e-4 |

The remaining temperature difference is a separate conversion issue. Pinned
`PeleCUtilities.H::pc_ctoprim` derives primitive energy from total energy minus
kinetic energy; PeleF's held MOL primitive conversion uses the stored internal
energy. Spatial interpolation generally makes these different even when both
interpolated conserved quantities agree. Changing this choice needs to retain
the existing dual-energy, source-funding and rollback contracts.

## Reproduction and evidence

Artifacts: `build/eb-fine-halo-diagnosis-20261019/`.

- `before/`, `source/`, `increment.patch`, `before-identity.json` and
  `source-identity.json`: this increment in the pre-existing dirty `main`
  checkout, based on `0a0e00ebe9138145fbb0e0e745665d0d768df009`.
- `diagnostic-build-commands.json`, `diagnostic-run-commands.json`, native
  `spatial_*.csv`, `fortran_spatial.csv`, `linear_spatial.csv` and
  `linear-halo-comparison.json`: diagnostic sources, exact commands and values.
- `check-commands.json` and `check-results.json`: focused regression,
  fixed-reference trajectory and reacting application commands/results.
- GNU Fortran 13.3.0, Debug checks and floating-point traps, Open MPI 4.1.6,
  CMake 3.28.3; existing `build/reacting-hot-diagnosis-20261019/registration`.
  MPI uses `/usr/bin/mpiexec`, with one OpenMP thread per process.

The first build attempt exposed an inappropriate AMR module dependency in the
lower transfer library. The planner now accepts integer layout metadata, so
no build-system dependency change is needed. The subsequent implementation
build completed; the final build includes input-shape admission checks and
focused affine/blocked-stencil/periodic/limiter/rollback assertions.

Final verification completed with the recorded source and binary hashes:

- Mechanism provenance/regeneration, Python syntax, `git diff --check` and
  `check_project_contract.py --source .`: pass.
- Focused CTests: 2 passed, 0 failed, 0 skipped. Transfer assertions take
  0.07 s; three-rank native MOL partial recovery takes 13.19 s.
- Hot CVODE subcycling application: pass (38.18 s), including 5 successful
  runs and 2 expected invalid-input refusals. Serial, MPI2 and MPI2-to-MPI3
  restart gas snapshots agree byte-for-byte, including intermediate outputs.
  Evidence: `application-runs/reacting_restart_n23__nhu/checks.json`.
- Fixed-reference comparison: all 3 steps complete (42.89 s), but 20 field
  limits fail; the checker returns 1. Evidence:
  `reference-runs/matched_eb_amr_6mu9dg03/checks.json`.

At 300 ns, maximum density difference is `1.871561e-7 kg/m3`, x-momentum
`6.557755e-4 kg/m2/s`, and temperature `2.781487e-3 K`. Before this increment
these were `1.811904e-7`, `6.559828e-4`, and `2.778743e-3`, respectively.
Thus the initial halo correction does **not** materially reduce the endpoint
maxima; density and temperature are slightly worse. The temperature limit
remains `1e-3 K`. There is no change to any acceptance limit or reference hash.

Maximum mass/element relative drift is `5.968818e-14`; energy-budget residual is
`2.413509e-17`. Hydrogen-consumption relative difference is `2.241819e-5`.
These conservation results do not substitute for trajectory agreement.
The next comparison must isolate the primitive energy choice and fine-level
stage evolution while preserving the paired gas/held-energy contract.
