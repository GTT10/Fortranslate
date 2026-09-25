# Serial EB tree periodic halos — 2026-09-26

Source `0145424`, parent `f02cb97`, local main. The serial halo constructor
accepts paired periodic faces and checks opposing root face metrics. Integer
donor coordinates wrap before peer/ancestor lookup for both face-only and wide
halos. Nonperiodic coordinates retain existing inflow/outflow handling.
Periodic hydro, transport and recursive updates explicitly reject the mode
until periodic seam flux reconciliation is implemented.

## Completed checks

Existing GNU Fortran 13.3.0/Ninja Debug tests-enabled build, Python 3.12.3.
No compiler or backend changes.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_halo_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_halo_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 passed**, 0.21 s test / 0.26 s total, no failures or skips. Added checks
cover a boundary child's opposite-side peer donor, six-layer root halos with
multiple-domain corner wrapping, and unchanged output arrays on unpaired
periodic labels or incompatible root EB faces. Existing prescribed inlet,
Mach inlet, restriction and temporal interpolation checks also pass.
Provenance, Python syntax, project contract and whitespace checks passed.

Artifacts: `build/eb-tree-periodic-halo-checks/build.log` and `LastTest.log`.
Mechanism SHA:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

This qualifies serial donor mapping only. Fine-level periodic seam metric
compatibility and conservative time updates remain future work. MPI periodic
tree exchange, application input, full configuration matrix, installed audit,
external solver parity and physical validation were not tested here.
