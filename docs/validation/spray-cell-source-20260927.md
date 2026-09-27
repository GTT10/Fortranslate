# Common particle SDC cell deposition — 2026-09-27

Implementation: `63f80152`.

Pinned PelePhysics `SprayParticles.cpp` interpolates gas values from neighboring
cells but deposits particle sources into `ijkc`, the cell containing the particle.
The recursive Fortran SDC route already selects this single-cell rule. The
common-step route still used weighted deposition to neighboring cells.

The existing `containing_cell` selector now passes through the gas deposition
and mixture-coupling wrappers; common particle SDC explicitly enables it.
Ordinary coupling retains its default weighted deposition. Sampling, particle
updates, integrated source magnitudes, cut-cell volume normalization and
transactional publication remain in their existing kernels.

Common-SDC spray checkpoint identity now includes `COMMON_SDC_CELL_SOURCE_V1`.
Snapshots made with the previous deposition rule are rejected rather than
silently continuing with changed numerics.

## Completed checks

Source-pinned fitted heptane provenance passed. GNU Fortran 13.3.0 Debug,
OpenMPI 4.1.6, MPI/spray ON, tests/SUNDIALS OFF, existing
`build/heptane-simple-mpi-debug` configuration:

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
python3 tools/prepare_evaporation_reference.py --bundle build/heptane-3sp-simple-reference/heptane_3sp.json --output build/spray-evaporation-cell-source
# From the prepared directory:
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none /usr/bin/mpiexec -n 1 /home/gtt13/projects/pelec/build/heptane-simple-mpi-debug/pelef_mpi_reactive_eb_tree_3d_selected_spray tree.nml spray.nml
```

Build completed; log `build/spray-cell-source-build.log`. Three application runs
completed: two-step evolution, one-step checkpoint, and one-to-two-rank restart.
Restarted gas CSV is byte-identical to the continuous run; particle rows agree
exactly after combining rank outputs by ID. One old-deposition restart was
rejected with `EB spray context mismatch` and no result CSV. Artifacts:
`build/spray-cell-{checkpoint,restart}/` and
`build/spray-evaporation-cell-source/restart-checks.json`.
Composite absolute inventory residual: `2.7105054312137611e-20`.

## Reference intervention

The unchanged pinned C++ two-step evaporation reference and its source identity
are documented in the [full reference record](evaporation-full-reference-20260927.md).
The immediately preceding Fortran run used the same Euler update and Simple
transport, with weighted deposition. Only the deposition route changes here.

| Gas field | Previous normalized maximum difference | Current difference |
| --- | ---: | ---: |
| Fuel density | `2.000e-4` | `1.395e-11` |
| Total density | `1.130e-9` | `6.060e-13` |
| Axial momentum density | `1.899e-4` | `5.156e-6` |
| Total energy density | `1.890e-7` | `2.372e-10` |
| Stored temperature | `2.270e-9` | `2.218e-9` |

Each difference is normalized by the maximum absolute reference value of that
field over all 64 cells. Particle temperature differs by `3.979e-13 K`, axial
velocity by `3.930e-14 m/s`, and mass by `2.019e-27 kg`.

This intervention identifies the deposition stencil as the dominant preceding
fuel-density discrepancy. Momentum and stored-temperature differences remain;
their cause is not established. Full gas/particle parity is not claimed.

```sh
python3 tools/compare_evaporation_particles.py --reference build/spray-evaporation-sdc-run2/spray00002.p3d --candidate build/spray-evaporation-cell-source/fortran.csv.particles.rank0.csv --reference-gas build/spray-evaporation-euler/reference-gas.csv --candidate-gas build/spray-evaporation-cell-source/fortran.csv --output build/spray-evaporation-cell-source/full-comparison.json
```

The diagnostic records input hashes and has no acceptance gate. Python syntax,
project contract and diff checks passed. No new CTest matrix, installed-binary
audit or physical validation was run.
