# Serial common-step EB tree viscosity — 2026-09-25

## Source and scope

Code `b412a28`, parent `83f057c`. The serial EB tree halo API optionally returns
volume fractions and completes interior, edge and corner donors. Existing
open-face donors and peer aperture/centroid checks remain. Donor search uses
same-level nodes first, then ancestors; physical coordinates clamp to the
domain. Ancestor state/temperature interpolation uses supplied time fractions.
Covered donors stay zero without EOS evaluation. Multiple active physical
inlet corners select the first axis in x/y/z order. Outputs publish together.

Common-step serial EB tree hydro accepts a finite nonnegative explicit
coefficient and requests this halo for positive values. It forwards coefficient
and volume fractions into each Euler step before existing peer reconciliation,
restriction, reflux and boundary/wall accounting. Omitted/zero retains the old
face-only halo path. Recursive/coupled and MPI wiring remain outstanding.

**Geometry limitation:** ancestor volume fractions use PCM donor values;
this does not reconstruct true fine-scale ghost geometry. Conservation checks
do not establish agreement with PeleC's geometry halos or its multi-inlet
corner policy. Reference geometry fidelity remains required work. This is
intermediate integration evidence, not complete EB-port or physical-validation
evidence.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug, Python 3.12.3; reused
`build/spray-amr-reactive-debug`. Focused CTest **2/2 passed**: existing zero
coefficient 1.06 s; coefficient 0.037 variant 1.14 s; total 2.25 s. No failures
or skips. Compiler warnings remain.

The existing fixture covers 4/5-node nested peer trees, axis planes in x/y/z,
outflow and six-face Mach input (including prescribed pressure), nontrivial
deep updates, inventory balance including root boundary/exposed wall terms
(existing scaled 3e-14 bound), and invalid-control/solver/step rollback.
No tolerance changes. This is bounded hierarchy evidence; it does not separately
prove temporal interpolation of positive-viscosity corner halos.

Mechanism provenance, Python syntax, project contract and whitespace checks
passed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_hydro_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_(hydro|viscosity)_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/eb-tree-viscosity-checks/build.log` and copied `LastTest.log`.
