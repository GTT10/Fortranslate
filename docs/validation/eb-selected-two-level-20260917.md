# Selected two-level EB application — 2026-09-17

Implementation: `e6bcab58e21a2cfe27cc81de22d6e943e8cf6fa4`, based on
`db795e9` on local `main`.

`pelef_reactive_eb_two_level_3d_selected` uses the same application source as
the fixed frontend. CMake configures its bundle imports and enables the
selected loading branch; the target links the existing selected EB runtime
and participates in installation. No generated mechanism rates were edited.

Named input mole fractions resolve against the bundle. The existing runtime
validates thermo/reaction/transport consistency and the common temperature
range for both initial density states. The default reaction integrator comes
from the bundle; the input may select explicit or implicit integration.
Selected checkpoint patch records use the existing schema with bundle SHA,
resolved mole fractions and integrator identity. The composite reader still
validates every patch and budget before publication. Unallocated optional
context retains the original fixed-checkpoint schema in the fixed executable.

Both one-patch common/subcycled stepping and multiple-patch common stepping
are connected. Element inventories currently recognize only H2, H, O, O2, OH,
H2O, HO2, H2O2, N2 and AR, and other species are rejected even with chemistry
disabled. This does not establish carbon-fuel EB support, CVODE integration,
multiple-patch subcycling, dynamic/multilevel EB or MPI EB.

## Independently executed checks

GNU Fortran 13.3.0, Debug/Ninja, existing `build/selected-tree-mpi-debug`,
configured with `/usr/bin/mpifort` (OpenMPI 4.1.6). These regressions are serial;
compiler/MPI/backend settings were not changed.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug \
  --target pelef_reactive_eb_two_level_3d pelef_reactive_eb_two_level_3d_selected --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug \
  --target test_selected_full_reactive_eb_two_level_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug \
  -R '^(regression_selected_eb_two_level_3d|regression_reactive_eb_two_level_3d)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tools/check_selected_eb_two_level_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, zero failures/skips: fixed regression 7.62 s, selected 6.50 s,
14.15 s total. Provenance, Python syntax, contract and whitespace passed.
Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

The selected fixture uses the pinned ten-species H2/O2 bundle, implicit
chemistry, transport and automatic interval selection. Two adjacent fine
patches and a single subcycled fine patch each complete three steps to 3e-8 s
with no rejected trials. Normalized composite residuals are respectively
**2.5359209340749391e-16** and **5.5147128964326066e-16**, with minimum limiter 1.
Each output has 104 active cells and ten species. Independent CSV integration
matches final inventories; values are finite, density/temperature positive,
species close mass, and chemistry has nonzero species contributions.

For each topology, one-step checkpoint continuation produces state and budget
CSVs byte-identical to uninterrupted execution. Corrupt bundle SHA and changed
composition reject restart without overwriting existing outputs. An unknown
input species also rejects initialization. The fixed regression retains its
earlier retry, subcycling, multiple-patch and corruption checks.

Logs: `build/eb-selected-two-level-checks/{build,test-build,ctest}.log`.
Selected inputs/checkpoints/CSVs:
`build/selected-tree-mpi-debug/tests/selected_eb_two_level_3d/`.
No full configuration matrix, clean installation audit, external solver parity,
convergence study or physical validation was run for this change.
