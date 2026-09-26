# Recursive physical-particle source synchronization — 2026-09-27

Implementation `86972b3544df861c46dd572816f1e7e9087d50f2`, based on `f7fa832e`,
local main, version 0.246.0. This resolves the preceding two-level
[crossing regression](spray-recursive-crossing-20260927.md).

## Numerical scope and reference boundary

Pinned PelePhysics `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b` deposits sources
into the particle-containing cell and transfers same-level source halos. PeleF
keeps that provisional physical/virtual/ghost stage source and its hydro
predictors. It additionally records the actual reaction from each physical
half kick, at its pre-drift position, on the finest containing physical cell.
The virtual and ghost copies are excluded from this independent field.

At the outer synchronization, the cellwise correction is:

```
actual physical half-kick increments on composite leaves
  - accepted provisional particle source increments on those leaves
```

Each old and final new half kick is counted once. Intermediate SDC iterations
return the cached old physical rate and are not counted again. Provisional
source increments are removed only at their own logical level. The correction
is applied to the private gas fields on uncovered fluid cells; temperature is
recovered with the existing EOS checks and parents are averaged down. The
particle-exchange ledger receives the same integrated correction. Wall/outflow,
chemical and external ledgers are not reclassified. No uniform/global offset,
source multiplier, tolerance change or clipping is used.

This synchronization is a **PeleF conservation extension**, not a claim that the
pinned upstream performs this correction. It changes accepted source placement
and does not establish upstream parity or temporal order. It is applied after
the recursive gas step, so strongly coupled cases can fail EOS admissibility;
those failures retain all caller state. Further accuracy and three-level coupled
qualification remain open before recursive application integration.

## Implementation boundaries

- The composite deposition API accepts optional `containing_cell`, retains
  owner exchange, fluid-volume normalization, MPI mode agreement, and atomic
  publication. Existing callers retain their default weighted stencil.
- `physical_source` exports are optional through the frozen kick, combined
  population stage, SDC stage and recursive coordinator. Frozen/combined stages
  export increments per fluid volume; SDC/recursive stages export rates.
- SDC work caches old physical sources and agrees export mode across stages.
  Failed physical export or a later virtual/ghost failure publishes neither the
  particle populations nor the exported source.
- Source accumulation checks finite inputs and arithmetic bounds. The final gas,
  temperature, physical populations, reaction fields and output ledgers remain
  within the existing outer transaction. Existing reaction-source history stays
  the gas SDC solve's result; no post-sync chemistry solve is inserted.

## Independently completed checks

Artifacts: `build/spray-recursive-source-sync-checks/`.
Existing Ninja Debug fixed-mechanism build `build/spray-amr-reactive-debug`,
GNU Fortran 13.3.0, `/usr/bin/mpifort` and `/usr/bin/mpiexec` OpenMPI 4.1.6.
Builds succeeded with existing real-comparison and executable-stack warnings.
Provenance, Python syntax, contract and whitespace checks passed. No failed
builds, failed tests or skipped tests occurred in this change's checks.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_reactive test_spray_mpi_eb_sample --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_mpi_eb_sample_2|spray_mpi_eb_sdc_recursive(_crossing|_entering)?_2)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **4/4 passed**, 32.55 s total:

| Case | Time | Maximum relative composite residual |
| --- | ---: | ---: |
| physical source export / population stages | 0.75 s | checked separately against physical kick |
| interior recursive coupled case | 17.28 s | 1.1112755134928231e-15 |
| fine to coarse crossing | 7.24 s | 1.5174158292877634e-15 |
| coarse to fine crossing | 7.23 s | 1.5176124431936805e-15 |

The export test has physical/virtual/ghost copies at the same location and
checks that only one physical reaction reaches composite leaves, including a
fine ghost position whose receiving cell is coarse. A third-population failure
retains the previous physical export and all population/source outputs.
Both crossing tests check final physical ownership, unchanged ID, heating,
evaporation, six half kicks and [1,2] level counts, with chemistry and all
molecular transport flags. The unchanged composite tolerance is **3e-12**.
Existing hydro rejection and post-level callback rejection exercise outer
rollback; final EOS-correction rejection was not independently injected.

Before arithmetic guards and reverse-crossing checks were added, the forward
crossing passed once (7.16 s, `ctest-crossing.log`); this is not counted again
as independent final evidence. `build-export.log`, `build-sync.log`,
`build-final.log`, `ctest-final.log` and `LastTest-final.log` preserve the run
history. Full matrix, changed-rank restart, three-level coupled evolution,
physical validation and external solver parity were not run.
