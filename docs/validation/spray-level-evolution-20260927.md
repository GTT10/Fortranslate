# Particle evolution in a capped AMR gas view — 2026-09-27

Implementation `6ea42a54`, based on `669aba22`, local `main`, version 0.246.0.
Pinned PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Implementation

Optional `maximum_level` in `advance_mpi_amr_mixture_coupling` reaches the
initial gas validation/deposition, gas sampling, spatial substep bound,
accepted gas deposition and final regular/EB parcel migration. The gas/EOS
deposition wrapper now forwards the ceiling to raw deposition, checks the
candidate gas and restores finer temperature arrays exactly before publication.
Gas and feedback above the ceiling remain unchanged and zero, respectively.

An explicit ceiling selects provisional level-view semantics: parent
restriction is skipped even at the finest ceiling. Omission preserves the
original composite update. Capped calls with wall-film populations reject;
their level-aware ownership remains unconnected. This API can evolve coarse
virtual copies, but does not select populations or control recursive clocks.
Fine ghost evolution still needs the halo operators wired into the stepper.
No recursive application or PeleC scheduling parity is claimed.

## Completed verification

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6. `spray_mpi_amr_step_2` and `_4` passed
(2/2, 0.86 s and 0.78 s; total 1.69 s).

The added fixture caps a three-level tree at level zero, with deliberately
different fine gas velocities. It compares against a single-level root tree
using the same initial particles and gas. Parcel positions, velocity,
temperature, mass, components, IDs, root gas and temperature match exactly;
substep counts also agree. Finer gas and temperatures are unchanged exactly,
finer source arrays are zero, and exported sources reconstruct the gas change
within the existing 3e-13 scaled tolerance. Rank-disagreeing ceilings reject
and retain gas arrays. Existing composite conservation, late rollback,
periodic movement, heating/evaporation, drag and KH-child assertions also run.
The new capped path was exercised for regular geometry at level zero; capped
EB evolution and intermediate ceilings were not independently run here.

The initial build failed on an overlong Fortran argument declaration. Splitting
the line fixed it; `build-initial-failed.log` retains that failure. The final
target build, provenance, Python syntax, project contract and whitespace checks
passed. No test failures or skips occurred. No full matrix, installed audit,
external solver parity or physical validation was performed.

## Reproduction

Artifacts: `build/spray-level-evolution-checks/` including copied CTest details.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_step --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_step_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
