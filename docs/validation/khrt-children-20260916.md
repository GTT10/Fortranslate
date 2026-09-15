# KH shed-mass children — 2026-09-16

## Source and scope

Tested parent `3be79c8663db4ef3e7d249a04febcbc0a281ab6d` plus this change.
Exact source/test hashes and logs: `build/khrt-children-checks/`.
The reference manifest pins ReitzKHRT.H, SpraySB.cpp and SprayParticles.cpp
at PelePhysics `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`; their local pinned
blobs were SHA256-checked. Existing default representative breakup is retained.

Optional KH shedding accumulates N*r³, releases a distinct child radius, and
subtracts emitted mass from the parent. Four symmetric transverse children
preserve mean momentum; a gas-energy debit accounts for their kinetic energy.
CFD substeps, parent lineage, CSV state and checkpoint persistence are connected.
The reference normalization is explicitly supplied in SI units. See the
[spray guide](../spray.md#kh-shed-mass-children) for formulas and differences.

## Focused verification

GNU Fortran 13.3.0, Debug with bounds/FPE checks, MPI disabled:

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d test_spray_core test_spray_integration --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(khrt|integration_(rollback|checkpoint)|application_cli)$' --output-on-failure
ctest --test-dir build/port-spray-debug -R '^spray_integration_rollback$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The first four-entry run passed the KH kernel, checkpoint and CLI entries;
the new rollback test failed because its chemistry-enabled argument was false.
Correcting that test argument enabled the intended late invalid-backend failure.
The rollback entry then passed independently in 0.11 s. All four targeted
entries therefore have completed passing results; there was no full rerun.
CLI completed all 33 Python tests without skips in the 43.64 s CTest entry.
Provenance, Python syntax, project contract and whitespace checks passed.
Earlier development compilation needed missing breakup helper imports corrected.

Kernel checks cover accumulation, partial/full parent release, four-child mass,
momentum/energy closure and invalid-input rollback. Application checks cover
active release and accumulation-only cases, conservation and byte-identical
final gas/parcel CSVs and checkpoints after stop/restart. Changed reference and
incompatible model/splitting settings reject. The integration test first proves
children are created, then checks gas, parents, child count and shedding-state
rollback on a later chemistry failure.

## Remaining boundary

Four deterministic transverse directions differ from the upstream random
azimuth/count algorithm. The minimum radius uses the configured diameter floor.
Wall impact and computational splitting cannot yet be combined with this mode.
RT representative-size changes do not themselves release children. General
turbulent dispersion, collisions, mixture droplets and particle MPI/AMR/EB
remain open. These checks establish bounded implementation behavior, not
external trajectory parity or physical validation. No full matrix was run.
