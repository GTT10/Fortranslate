# Serial EB tree CTU subcycling and coupled composition — 2026-09-25

## Source and scope

Implementation: `70388ef`, based on `74cfda454fab9b371113e4cd43413d2d2bef3a94`.
PeleC reference remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
The optional extended geometry cache now selects PLM CTU in recursive serial
hydro and in both branches of the full reaction/transport composition.
Six-layer state preparation uses time-interpolated ancestor donors. Existing
peer reconciliation, fine-substep flux accumulation, reflux and restriction
handle the accepted integrated CTU flux. Covered storage is preserved.
Invalid cache/solver combinations reject without publishing private stages.

## Focused evidence

GNU Fortran 13.3.0, Python 3.12.3, Ninja Debug, existing fixed-mechanism
`build/spray-amr-reactive-debug` configuration. These are serial tests; no MPI
execution or external solver comparison was performed.

The five-node, three-level z-plane fixture has discontinuous densities and
1000 K initial fluid temperatures for CTU, within the mechanism EOS range.
Both explicit and implicit chemistry run all four chemistry/transport switches.
The fully coupled explicit result equals independently invoked
R/2, T/2, H, T/2, R/2 operators exactly, including ledgers. Deliberate hydro
rejection after private reaction/transport stages preserves state and diagnostics.
CTU uses `difmag=.037`; its initial velocity is zero, so this does not isolate
a nonzero artificial-viscosity correction. No tolerance or mechanism changes.

| CTest | Result | Seconds |
| --- | --- | --- |
| `unit_eb_tree_subcycle_full_3d` | PASS | 47.21 |
| `unit_eb_tree_ctu_subcycle_full_3d` | PASS | 143.14 |
| `unit_eb_tree_ctu_common_full_3d` | PASS | 95.53 |

3/3 passed, total 285.93 seconds; no skipped or interrupted runs. Maximum
scaled composite residual: recursive CTU `5.627e-15`, common CTU `5.958e-15`,
against the unchanged `3e-14` threshold. Initial compilation failed on one
new overlong test import line; corrected before executing tests. Compiler
warnings remain; this is not a warning-clean build.

Commands from the repository root:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_subcycle_full_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_(ctu_(subcycle|common)_full|subcycle_full)_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, Python syntax, project contract and whitespace checks passed.
Pinned mechanism SHA:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Artifacts: `build/eb-tree-ctu-subcycle-checks/build.log` (initial failure),
`build-final.log`, and `LastTest.log` in that same directory.

MPI CTU stepping, tree application selection/restart, broader physical-boundary
halos, complete reference parity and physical validation remain outside this check.
