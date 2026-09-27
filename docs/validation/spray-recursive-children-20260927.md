# Recursive SDC KHRT shedding children — 2026-09-27

Implementation `9eac6d3eca1472050a0343dfafa4b95ae50f6bab`, based on `0925e7b5`,
local main, version 0.246.0.

## Change

The recursive population APIs now carry a private optional `last_id` through
old-time physical stages. Their halo kernel samples stage-time gas/film
properties, applies frozen-gas exchange, then calls the existing collective
KHRT release/child constructor. The transverse child kinetic energy is debited
from gas at the original parent source position. Parent and child component
inventories then undergo the existing drift and level redistribution.

Each enclosing layer publishes its ID watermark only after its local operation
succeeds; the outer wrapper publishes only after gas synchronization, physical
collection and optional computational splitting succeed. Intermediate iterations
do not generate children. Virtual/ghost populations and the new-time half kick
use diameter evolution without child generation. KHRT's reference-volume control
is cleared in those private no-child options. Computational `split_created`
continues to count only multiplicity splits, not KHRT shedding IDs.

Reflecting/EB contacts take precedence over release, using a preview of the same
segmented flight as the eventual drift. Outflow remains part of the ordinary
boundary inventory. The halo kernel now also defers level-source publication
until its optional physical export succeeds, preserving the complete source
transaction along with population and ID publication.

The application no longer rejects `khrt_children`; existing option validation
still requires KHRT and a positive reference diameter cubed. Films remain
unsupported in recursive SDC.

This reuses PeleF's common-step KHRT release/drift convention. It is not a claim
of exact parity with pinned PelePhysics
`6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`: upstream child construction occurs
after the tile's particle-update loop, and detailed child-placement/time-order
parity remains to be established. A subsequent
[placement audit](spray-child-placement-20260927.md) resolves the interior
position formula and confirms that copies and new-time stages create no children;
other reference differences remain. This check does not qualify the physical
breakup model or its empirical constants.

## Independently completed checks

Artifacts: `build/spray-recursive-children-checks/`.
Existing fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`,
GNU Fortran 13.3.0, `/usr/bin/mpifort` and `/usr/bin/mpiexec` OpenMPI 4.1.6.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d test_spray_mpi_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc_recursive_children_2$' --output-on-failure
# Add an explicit transverse-energy assertion and rebuild the focused regressions:
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d test_spray_mpi_eb_reactive test_spray_mpi_eb_sample --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_mpi_eb_sdc_recursive_children_2|spray_mpi_eb_sample_2)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final **2/2 tests passed**, population-stage regression 0.75 s and recursive
children 8.05 s, 8.85 s total. No test failures or skips. Initial child-only check
also passed, 8.05 s. Provenance, final builds, Python syntax, project contract and
whitespace checks passed. One initial build failed at an overlong argument list;
it was wrapped before any test ran. That log is retained in
`build-line-length-failed.log`; existing executable-stack linker warnings remain.

The two-rank child fixture has two AMR levels, gas chemistry and transport, and
one fast mixed-liquid parcel. Drag, heating and evaporation are disabled to
isolate the release energy. It checks rejected hydro after old-time particle
work, child-capacity rejection, unchanged caller inventories/IDs on failure,
and a successful retry. The accepted result has **9 live particles, watermark
9**, globally unique sorted IDs on the expected owner, no computational splits,
and gas exchange only in energy: **-9.6937861994497204e-11 J**. The coupled
inventory residual is **1.0551308123124431e-15**, below the unchanged 3e-12
criterion. The second regression exercises optional-ID omission, stage ordering,
acceleration and late-copy rollback on the existing no-child path.

Application wiring was compiled. No application runtime, wall-contact release
case, multi-owner simultaneous shedding, accelerated/evaporating child case,
rank-changing restart, selected mechanism, three-level evolution or physical
validation was run in this focused check.
