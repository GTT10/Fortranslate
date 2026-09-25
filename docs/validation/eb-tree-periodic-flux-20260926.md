# Serial EB tree periodic flux synchronization — 2026-09-26

Source `c9411f7`, parent `090f4e9`, local main. The synchronization API adds
optional `boundary_faces`. Paired periodic root metrics are validated before
private changes. Same-level raw periodic transfers use the high-domain face
as owner; a copied owner flux supports a patch spanning the full domain.
Fine transfers are restricted and mirrored into the opposite periodic coarse
face, including the root itself, before coarse exposed-cell correction.
All effective fluxes and corrected states publish together after averaging.
Serial common-step hydro now forwards periodic labels through this path for
PCM and extended-geometry PLM/CTU updates.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug, Python 3.12.3; existing tests-enabled build.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_flux_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_flux_3d$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_hydro_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_hydro_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Two focused tests passed independently, no skips:

- Flux synchronization: 1/1, 0.31 s test / 0.37 s total. Three plane normals,
  four/five-patch trees, ordinary and periodic cases. Periodic grandchildren
  touch the domain seam, with one-sided and two-sided refinement. Root open-face
  transfers agree; the largest periodic composite inventory residual was
  7.46e-15 against the existing 3e-14 threshold. Nonfinite and nonphysical
  corrections preserve all tree and flux outputs.
- Common hydro: 1/1, 26.32 s test / 26.36 s total. PCM and PLM/CTU periodic
  updates pass the existing 3e-14 composite conservation bound in all six
  plane/tree combinations, alongside the existing ordinary and Mach checks.

The first flux fixture incorrectly compared unused flux values on covered
faces; it now compares aperture-weighted transfers. Its nonphysical-flux
fixture also targeted a newly covered grandchild after moving refinement to
the seam; it now targets an exposed periodic root face. Failure logs are
retained. No numerical tolerance or production validation was relaxed.
Provenance, Python syntax, project contract and whitespace checks passed.

Artifacts under `build/eb-tree-periodic-flux-checks/`: `build.log`,
`build-hydro.log`, `flux-LastTest.log`, `hydro-LastTest.log`,
`initial-failure.log`, `rollback-fixture-failure.log`.
Mechanism SHA:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

Periodic recursive updates, molecular transport, MPI time updates and tree
application inputs remain unconnected. No full matrix, installed audit,
external solver parity or physical validation was run.
