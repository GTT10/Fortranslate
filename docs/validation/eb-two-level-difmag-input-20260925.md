# Two-level/multipatch EB viscosity input — 2026-09-25

Code `511c7d3`, parent `952aa5e`. Shared fixed/selected application source accepts
finite nonnegative `difmag` in `&eb_two_level`, forwards it through single-patch
and multipatch full composition, and copies it into patch checkpoint config.
Positive values therefore bind restart identity. Default zero retains the old
solver identity; -1 is invalid, with no automatic sentinel in this input.

## Completed checks

GNU Fortran 13.3.0, Python 3.12.3, Ninja Debug; reused
`build/spray-amr-reactive-debug`. Fixed executable built. **2/2 CTests passed**:
single-patch 4.21 s and multipatch 4.54 s (8.79 s total). No failures or skips.
Provenance, Python syntax, project contract and whitespace checks passed.

Both reacting/transport planar fixtures use fixed root intervals and recursive
stepping. Coefficient 0.037 gives exact continuous/resumed state CSV equality;
interrupted output differs from the final reference. Changed 0.05, negative
and NaN values reject before new CSV publication. Omitted/zero outputs match,
and an omitted-coefficient checkpoint resumes with explicit zero exactly.
Multipatch automatic timestep selection is disabled in this focused fixture.

Selected builds, regridding/adaptation with positive coefficients, budget CSV
parity and historical executable restart were not independently run. These
are bounded application checks, not external solver parity or physical
validation. Uniform reference defaults and EB geometry fidelity remain open.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_eb_two_level_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^reactive_eb_two_level_difmag_(reference|multipatch)$' --output-on-failure
python3 -m py_compile tools/check_eb_difmag_input.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/eb-two-level-difmag-input-checks/build.log` and `LastTest.log`.
The test log records unique fixture paths for each layout.
