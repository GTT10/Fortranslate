# KHRT child placement and stage audit — 2026-09-27

Regression `47891c33b0ea4748cd4ed04f13c44466b835617a`, based on `ab96f409`,
local main, version 0.246.0. No solver arithmetic changed in this check.

## Reference finding

The preceding recursive-child record left child placement unresolved because
pinned PelePhysics creates children after its tile particle-update loop.
Inspecting the complete path at
`6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b` resolves the interior placement formula:

- `Source/Spray/BreakupSplash/ReitzKHRT.H`, around line 152, saves the parent's
  current location before parent drift, together with its velocity and release
  properties.
- `Source/Spray/SpraySB.cpp`, around line 179, initializes each child at the saved
  release location plus `sub_dt` times that child's velocity.
- Thus later allocation does not imply that the child remains at the parent's
  final location or receives no displacement. PeleF's construction at the release
  location followed by one full drift gives the same interior position formula
  for a given child velocity and interval.

This finding removes the suspected interior placement-order mismatch. It does
not establish identical children: upstream chooses their number using its
multiplicity exponent and samples random azimuths; PeleF currently uses four
symmetric directions with a compensating gas-energy debit. Exact child count,
angular distribution, source treatment, wall encounters and detailed subcycling
parity remain outside this evidence.

## Completed checks

Artifacts: `build/spray-child-placement-checks/`. Existing fixed-mechanism Ninja
Debug build `build/spray-amr-reactive-debug`; GNU Fortran 13.3.0 and OpenMPI 4.1.6
(`/usr/bin/mpifort`, `/usr/bin/mpiexec`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_sample --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sample_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final **2/2 passed**: 2 ranks 0.74 s, 4 ranks 0.77 s, total 1.57 s; no skips.
Provenance, build, Python syntax, contract and whitespace checks passed.
The added old/intermediate/new-stage fixture checks:

- only the physical population creates four children and advances the ID
  watermark from 1 to 5; identical virtual and ghost copies create none;
- each child's interior midpoint position equals `birth_position + dt*velocity`
  exactly, using the pinned reference formula;
- physical publication remains deferred at the old stage;
- intermediate and final stages allocate no further IDs;
- the new-time half kick leaves every midpoint position exactly unchanged.

Initial test attempts failed for two fixture reasons, retained in the artifact
logs. Subtracting two O(1) positions before subtracting a small displacement lost
precision against a 1e-18 residual threshold. The final check instead compares
position with the directly evaluated reference expression using exact equality;
no tolerance was loosened. After that correction, the new-time kick correctly
rejected an outward-moving parcel starting at the outermost center of a one-cell
halo: its next interpolation needed another ghost cell. The final fixture starts
inside the fine patch. The diagnostic run confirmed stage 2 rejection with the
ID watermark still 5. Existing production recursive stages use three-cell halos.

No full upstream executable comparison, wall-contact case, application run or
physical validation was performed.
