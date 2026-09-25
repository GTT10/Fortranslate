# EB transport stage donor preparation — 2026-09-25

Implementation `8cb10ec`, based on `76780e3` on `main`.
PeleC pin remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

`advance_reactive_eb_transport_3d` now accepts an optional
`eb_transport_halo_provider_3d` procedure. It receives the current stage U/T and
geometry and fills a one-cell donor halo and six face-enable switches before
each SSPRK2 Euler stage. The second call receives the first stage's recovered
state. Callers without the procedure retain the existing closed transport faces.
An optional integrated-flux result uses one-half of each Euler stage's integrated
flux. State, temperature and flux publish only after final thermodynamic recovery;
a second-stage donor failure leaves all three caller-owned results unchanged.
Disabled transport or a zero interval succeeds with zero integrated exchange.

The provider owns boundary labels, prescribed values and donor preparation. This
interface does not itself implement periodic seam reconciliation or thermal-wall
flux policies. Passing it through the coupled gas driver and applications remains
pending; this is not a completed physical-boundary application path.

## Focused checks

Existing Ninja Debug directory `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, Python 3.12.3, configured OpenMPI wrappers. No toolchain/backend changes.
Two tests passed, **1.63 s total**, with no failures or skips:

- `unit_reactive_eb_transport_3d`: 0.98 s; existing closed-boundary and coupled checks.
- `unit_eb_transport_halo_3d`: 0.60 s; existing checks plus full-fluid and cut-cell
  SSPRK2 updates using the regular boundary filler with prescribed inflow and
  outflow. Checks observe both provider calls and changed second-stage state,
  compare independently assembled Euler stages and weighted fluxes, close every
  component inventory against aperture-weighted boundary transfer, and exercise
  a provider failure at stage two. The inventory criterion remains `2e-12` scaled
  by initial component inventory.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_transport_halo_3d test_reactive_eb_transport_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_(eb_transport_halo_3d|reactive_eb_transport_3d)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, Python syntax, project contract and whitespace checks passed.
Mechanism SHA unchanged:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Artifacts: `build/eb-transport-stage-halo-checks/{build.log,LastTest.log}`.
Compiler warnings remain. The test's host-associated internal donor procedure
also produces a GNU linker executable-stack warning; this is a test-executable
warning, not a newly built application qualification.
No full matrix, installed-binary audit, application run, external-solver parity or
physical validation was performed for this increment.
