# Independent MPI particle source export — 2026-09-27

Implementation `bac68877`, based on `68f31f43` on local `main`.
PeleC pin: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`;
PelePhysics pin: `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`.
H2O2 source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Capability

The generic and EOS-checked MPI AMR deposition APIs, and
`advance_mpi_amr_mixture_coupling`, now accept optional `gas_change`.
This owner-local field holds the time-integrated particle feedback per fluid
volume, in the same conserved-component order as the gas. Divide by the accepted
interval to obtain an interval-average source rate. Output replaces the previous
field; it does not accumulate across separate calls.

Deposition accumulates the same routed contributions into a separate zero field,
then applies ordinary or fluid-weighted EB parent restriction. It never subtracts
two large gas states to infer a small source. The EOS wrapper publishes gas,
temperature and source together. Globally empty deposition returns a zero source
without changing the gas. Presence of the optional source output must agree
across ranks and is checked collectively.

The particle interval sums only accepted substep increments. Rejected trials,
substep exhaustion and late migration failure leave the caller's source intact,
alongside the existing gas/particle rollback. The field includes gas exchanges
from drag, heating, evaporation and the KH child-energy debit. Wall and exterior
inventories remain separate; this field does not include standalone film
thermal/evaporation evolution. Generic regular and nodal-EB trees share this path.

## Reference and remaining integration

The pinned `Source/Particle.cpp` was inspected using `git show` in
`/mnt/sn850/Projects/gtt2/external/PeleC`. Its `particleMKD` and `particleMK`
transfer separate old/new spray sources after their kick/drift stages.
Independent source export is needed before connecting those source families to
PeleF's SDC gas stages. This increment retains the existing particle interval
algorithm. It does not implement that kick/drift schedule or enable the EB spray
application's SDC selector; those remain open.

## Completed focused checks

GNU Fortran 13.3.0, Ninja Debug, fixed-mechanism tests-enabled build
`build/spray-amr-reactive-debug`; `/usr/bin/mpifort` and `/usr/bin/mpiexec`,
OpenMPI/OpenRTE 4.1.6. Existing compiler warnings remain.

- Mechanism provenance and all three requested targets built successfully.
- CTest **3/3 passed**, 1.95 s total: `spray_mpi_amr_step_2` 0.79 s,
  `spray_mpi_eb_step_2` 0.59 s, `spray_mpi_eb_deposit_2` 0.52 s.
- Regular three-level mixed-particle coupling retains gas-plus-liquid
  conservation within 3e-12 relative, including heating, evaporation, drag,
  coarse/fine movement, periodic wrapping, KH children and migration.
- The nodal-EB interval retains fluid-volume gas/particle/wall conservation
  within 3e-12 relative, with reflection, heating, evaporation and ownership
  transfer. Both interval tests verify owner-local exported fields and gas
  reconstruction within 3e-13 times `max(1,abs(final gas component))`.
- Interval failure retains sentinel source values. Existing substep exhaustion,
  interval disagreement, migration capacity and gas/temperature/parcel rollback
  assertions remain active; globally empty regular populations pass.
- EB gas deposition rejects a nonphysical final state and retains the previous
  source. A subsequent `1e-25 * conserved_state` integrated contribution leaves
  the gas bit-for-bit unchanged because it is below gas-state resolution, while
  its source field retains the composite inventory with maximum relative error
  **1.1479437019748900e-16** (bound 3e-14). This directly distinguishes source
  accumulation from subtraction of rounded gas states.
- Python syntax, project contract and `git diff --check` passed. No failed or
  skipped tests were counted.

Artifacts: `build/spray-sdc-source-checks/`: `provenance.log`, `build.log`,
`ctest.log`, `detail.log`, `contract.log`, `toolchain.log` and
`reference-particle.cpp`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug \
  --target test_spray_mpi_amr_step test_spray_mpi_eb_step test_spray_mpi_eb_deposit --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 \
  OMPI_MCA_hwloc_base_binding_policy=none \
  ctest --test-dir build/spray-amr-reactive-debug \
  -R '^spray_mpi_(amr_step|eb_step|eb_deposit)_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

No new application/restart, selected-mechanism build, full configuration matrix,
installed-binary audit, SDC order study, reference trajectory parity, scaling or
physical spray validation was run. These are source-export and bounded coupled
conservation checks.
