# Regular 3D patch-tree conservation foundation — 2026-09-17

Parent: `b34534d`, local `main`. See the
[design and API boundary](../design_decisions/amr-patch-tree-3d.md).

Implemented multiple child boxes and arbitrary tree depth, checked topology,
generic conserved-state storage, PCM initialization, bottom-up restriction and
composite inventory. Existing two-level PCM/restriction kernels are reused.
The new module is part of `pelef_core`; existing applications are unchanged.
No tree flow evolution, ghost exchange, reflux, dynamic remap, MPI, EB or
application/restart integration is claimed by this check.

## Commands and configuration

Reused `build/spray-mpi-migration-gnu-debug`, GNU Fortran 13.3, Debug with
bounds/FPE checks (MPI and experimental spray configured). This change adds a
core source/test; compiler, MPI and backend configuration were unchanged.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target test_amr_patch_tree_3d test_amr_hierarchy_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^unit_amr_(patch_tree|hierarchy)_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Final focused tests: **2/2 passed, 0 failures, 0 skips**, 0.04 s.
The new test uses four levels/six nodes, adjacent root children and two refined
branches, level ratios 2/3/2, nonuniform two-component root fields and a deep
leaf perturbation. It checks physical-volume-weighted conservation, absence of
double counting, correct ancestor sampling, propagation to the root and the
unrefined special case. Invalid overlap, parent cycles, inconsistent level
ratios, wrong parent shapes and nonfinite data leave destinations unchanged.
The original two-level hierarchy unit also passed.

The first test exposed naive accumulation error of
`-1.1596057447604835e-11` in the first inventory component. Compensated
summation resolved this without changing the `1e-11` test threshold. The
initial failed attempt is in `ctest.log`; final results are in `ctest-final.log`.
Build, provenance, project contract, Python AST syntax for two verification
tools, and whitespace checks passed. No full matrix, installed-binary audit,
external solver parity, scaling or physical validation was run.

Logs and source hashes: `build/amr-patch-tree-3d-checks/`.

```text
cd9f847cc2e6341743b44b8222b201278c8189b564cb6caf0b8873539bbbb197  src/amr/amr_patch_tree_3d_mod.F90
c8286c1ce6788fbe37dd8497782840fa9403db8b9c59141ca36f28d55d052f76  tests/unit/test_amr_patch_tree_3d.F90
```
