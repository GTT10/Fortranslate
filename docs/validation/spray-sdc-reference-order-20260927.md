# Corrected SDC particle kick ordering — 2026-09-27

Implementation `a74881c4`, based on `9902c8ad`, local `main`, version 0.246.0.
Pinned PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Reference and correction

The pinned `Source/Particle.cpp::particleMKD` returns unless sub_iteration is
zero. `particleMK` returns unless sub_iteration is sub_ncycle-1 (apart from the
separate sub_ncycle=0 path). In `Source/Advance.cpp::do_sdc_iteration`, the first
iteration copies old sources into new sources before hydro. Thus intermediate
SDC iterations retain the old spray source; they do not reevaluate a new kick.
Both files were read using `git show` at the pinned commit and saved in the
artifact directory. The external checkout working tree was not the reference.

The previous level controller and common-step wrapper reevaluated a private new
half kick each iteration. Although this avoided repeated physical kicks, the
intermediate gas feedback differed from the pinned source scheduling. This
change corrects both `advance_mpi_eb_population_sdc_stage` and
`advance_mpi_eb_spray_sdc`: intermediate iterations return the saved old source,
while only the final iteration evaluates the new half kick. The level controller
also retains its old halo remainder and physical boundary ledgers. No new gas
sampling or particle update occurs during its intermediate callback.

The previous [level stage record](spray-population-sdc-20260927.md) describes
the historical repeated-evaluation implementation and its test result; its
intermediate-stage policy is superseded here. Final publication and rollback
remain transactional. This corrects scheduling, not the still-open recursive
population lifecycle, full source-halo assembly or application integration.

## Completed checks

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6. Three focused tests passed with no failures/skips:

- `spray_mpi_eb_sample_2`, `_4`: 0.68/0.67 s.
- `spray_mpi_eb_sdc_2`: 5.76 s.
- CTest total: 7.16 s.

The level fixture now verifies exact old-source/remainder retention at an
intermediate stage, even with an invalid cached ghost node that a new kick
would reject. That node does reject at the final stage, preserving the prior
rates and stage cursor. Restoring it permits retry; final rates equal the
independent new half-kick result divided by its actual duration, and particle
states match one old drift plus one final new kick.

The existing common-step reacting/transport/mixed-particle test passes its
coupled inventory, late splitting failure and rollback checks under the corrected
policy. This is a regression check, not external solver numerical parity.
Film-specific cases, the full configuration matrix, installed audit and physical
validation were not rerun. Provenance, Python syntax, project contract and
whitespace checks passed.

## Reproduction

Artifacts: `build/spray-sdc-reference-order-checks/`, including pinned reference
files, build, provenance, contract, CTest summary and copied `LastTest.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_sample test_spray_mpi_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_mpi_eb_sample_[24]|spray_mpi_eb_sdc_2)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
