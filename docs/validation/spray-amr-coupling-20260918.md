# Composite AMR spray field transfer — 2026-09-18

Implementation: `e373d15e502149e4ad505fbc973220065494063e`, based on
`68dc82e`. `spray_amr_tree_coupling_mod` builds a cloud-overlap stencil, samples
generic tree fields and deposits integrated contributions into volume densities.
The cloud width is the finest containing cell width in each axis. Only exposed
cells contribute; covered parents are excluded. Periodic cloud pieces wrap and
nonperiodic clouds clip/renormalize. Positive weights sum to one. The same
weights drive sampling and deposition; deposition divides by receiving volume
and refreshes covered parents with the existing bottom-up averaging kernel.

On uniform meshes this is ordinary trilinear cloud-in-cell transfer. Across
coarse/fine interfaces it treats fields as piecewise constant inside each
exposed cell; it is not a linear-exact reconstruction. This bounded conservative
transfer is not a claim of pinned PeleC AMR-particle parity or convergence.

Stencil arrays, samples and deposited state publish only after success. A
configurable stencil capacity (up to one million entries), invalid positions,
metadata and unsafe density arithmetic reject without committing. Empty
populations still validate geometry; an empty deposit leaves all fields unchanged.
This is a generic field API, without EOS/positivity recovery, owner-local MPI
field exchange, coupled trajectory stepping or EB solid treatment.

## Completed focused verification

GNU Fortran 13.3.0 Debug/Ninja, `/usr/bin/gfortran`. Reused the compatible serial
fixed and selected H2/O2 build directories established for AMR partitioning;
compiler, MPI and optional-backend settings were unchanged. CMake regenerated
the new source/test targets.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-partition-debug --target test_spray_amr_tree_coupling --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-partition-debug -R '^spray_amr_tree_coupling$' --output-on-failure --output-log build/spray-amr-coupling-checks/ctest-final.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-partition-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final **1/1 passed**, no failures/skips, test time **0.00 s**, total **0.03 s**.
An initial run also passed; empty-population validation and density-overflow
cases were then added and the focused test rerun. This is one unique test.

The fixture compares uniform sampling/deposition against the pre-existing
eight-point periodic particle stencil. A three-level tree checks constant-field
preservation, independently calculated piecewise-field samples of
`17, 14, 27.29, 14.64` across coarse/fine, nested and periodic boundaries,
positive normalized weights, receiving volumes and exclusion of covered cells.
Composite volume integrals after deposition equal the initial inventory plus
the supplied signed contributions within 3e-14 of the inventory scale. The
test also checks stencil-capacity rollback, invalid-position sample/deposit
rollback, clipped physical corners, empty populations and rejection of density
overflow before division.

The selected library compiled successfully. Provenance, Python syntax, project
contract and whitespace checks passed. No build/test attempt failed. No Python
source changed. This is serial field-transfer unit and selected-library compile
evidence, not a selected executable, MPI transfer, install audit, AMR convergence
study or physical spray validation.

Logs: `build/spray-amr-coupling-checks/{build,rebuild,ctest,ctest-final,build-selected}.log`.
Executable: `build/spray-amr-partition-debug/tests/test_spray_amr_tree_coupling`.
