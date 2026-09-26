# Recursive spray cell source and unresolved interface budget — 2026-09-27

Implementation `d7854311`, based on `2d087ccd`, local main, version 0.246.0.
**The new crossing regression fails. This change is not a conservation fix.**

## Source discrepancy corrected

Pinned PelePhysics `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`,
`Source/Spray/SprayParticles.cpp` samples gas with interpolation weights but
sets `cur_indx=ijkc` for source deposition. `ijkc` is the containing cell from
`floor((position-plo)/dx)`. Source normalization uses its cell volume, with a
volume-fraction divisor for a cut cell. Source addition occurs before drift.
`SprayParticles.H::transferSource` sums same-level overlaps and copies source
halos in its PeleC branch; it does not directly route unmatched ghost sources
into coarser valid cells.

The recursive physical/virtual/ghost population stage previously reused CIC
sampling weights for deposition. It now selects containing-cell deposition,
with finite/positive fluid-volume checks, overflow protection, MPI-agreed mode
and transactional output. Gas interpolation remains unchanged. The generic
halo API retains its existing weighted default for other callers; this change
does not establish reference agreement of the common-step coupling path.
Pinned source copies are saved in the artifact directory below.

## Unresolved crossing observation

The new two-level, two-rank planar EB regression starts a mixed parcel at
x=0.00199999 m, moving right across the fine-patch boundary x=0.002 m during
1e-8 s. Chemistry and all molecular transport flags are enabled. It verifies
accepted level counts [1,2], six half kicks, unchanged particle ID and final
root ownership before checking composite gas/particle conservation.

- Original weighted deposition: relative composite residual
  **2.2141349802987069e-9**, above the unchanged **3e-12** criterion.
- After reference cell deposition: residual **5.8168860550643627e-6**;
  the failure remains and is larger.
- With cell deposition, maximum absolute gas inventory residual after subtracting
  gas boundary/wall/chemical/external and deposited particle ledgers is
  **1.5751057418948996e-17**. This maximum mixes conserved-variable units and is
  diagnostic only, not a normalized acceptance metric.
- Particle change plus deposited exchange has energy mismatch
  **4.5875029072855817e-9** and mass mismatch **-5.7901786274070718e-16** in
  the fixture's SI inventory units.

This localizes the observed discrepancy to physical-particle versus composite
source accounting. Different coarse-copy and fine-physical time histories and
valid/ghost support are candidate causes; their individual contributions have
not been isolated. It does not prove that upstream has the same failure.
No global correction, source rescaling, relaxed tolerance, disabled test or
expected-failure override was used. The failing CTest remains registered.
Next work must resolve interface exchange/synchronization and preserve reference
provenance before enabling recursive application stepping.

## Completed checks and failures

Artifacts: `build/spray-recursive-crossing-checks/`.
Existing fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`,
GNU Fortran 13.3.0, `/usr/bin/mpifort` and `/usr/bin/mpiexec` OpenMPI 4.1.6.
Builds succeeded (existing compiler/linker warnings remain). Provenance preflight,
Python syntax, project contract and whitespace checks passed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_sampling test_spray_mpi_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_mpi_amr_sampling_2|spray_mpi_eb_sdc_recursive_2)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final focused check: **2/2 passed**, 17.30 s total. Sampling test 0.49 s checks
off-center contributions in exactly one fine ghost cell, periodic image mapping,
fluid-volume scaling and solid-containing-cell rejection despite neighboring
fluid donors, retaining the previous output. Interior recursive core 16.76 s
retains its original conservation and rollback checks.

Other runs (not counted again as independent final passes):

- Initial crossing check: 0/1 passed, 6.95 s; `ctest.log`, `LastTest-first.log`.
- Diagnostic crossing repeat: 0/1 passed, 7.21 s; `ctest-diagnostic.log`.
- After cell-source change: population stage test passed 0.75 s; crossing test
  failed 6.94 s; **1/2 passed**, 7.75 s total, `ctest-cell.log`, `LastTest-cell.log`.
- Final logs: `ctest-final.log`, `LastTest-final.log`; no skips in these runs.

The failed crossing outputs include OpenMPI/X authorization noise, but both ranks
reach the explicit composite-inventory assertion and MPI abort. This is a
numerical regression failure, not an approval or launch block.
