# Native held MOL root-window replay, 2026-10-11

## Implemented scope

The serial and MPI recursive gas APIs accept `native_replay_history` for a
complete ordinary (`redistribution="none"`) chemistry-enabled MOL root window.
They restore the complete gas/T/UEINT, reaction/heat and numerical-prefix clock
anchor and the private funding context. Native external force/work and auxiliary
stages are frozen; hydro, molecular transport, chemistry, face families and
funding receipts are rebuilt through the original recursive driver. Classified
spatial source providers retain their accepted cache rather than repeating their
callbacks. All output fields, fluxes, history and budgets publish only after the
whole window accepts.

Whole-window preflight checks the physical level clocks, refinement intervals,
ancestor time fractions, owner-local payloads, policies and source-cache shapes
before invoking the physics. Prepared redistribution, routed sources, sub-parent
windows, incomplete source prefixes and weighted native replay are outside this
entry. Automatic weighted parent repair and the actual future-fine producer
remain open. This library feature does not establish pinned PeleC trajectory
parity or physical validation.

## Source and artifacts

Work continues in the existing dirty `main` checkout at
`0a0e00ebe9138145fbb0e0e745665d0d768df009`; earlier user changes are preserved.
The branch contains the stored `origin/main` base
`5a2c818899f7213dad2fbb3c6808558711aced10`.
Source snapshots, incremental diffs, build logs, linked identities, CTest XML and
resume metadata are in `build/native-mol-parent-replay-20261011`.
The toolchain is GNU Fortran 13.3.0 / OpenMPI 4.1.6, Debug, with the fixed and
selected `h2o2_full` mechanism, spray enabled and SUNDIALS disabled.

## Verification status

The final candidate passes four CTest cases: the new native root replay test
(106.15 seconds) and the existing three-rank MOL, MOL reflux and held reset
controls (16.51 seconds together). A separately linked, unchanged external-stage
regression passes both one/two correction cases (25.48 seconds). The mechanism
provenance unit suite passes 14/14, 153 Python files pass syntax parsing, and
mechanism provenance, project contract and `git diff --check` pass.

The new test uses the fixed full H2/O2 mechanism with pure N2 and native forced
chemistry, cut geometry, nonzero body forcing and auxiliary gas energy. It checks
one/two corrections, seven `[1,2,4]` clocks, distributed `[0,2,0]` and owner-zero
layouts, full owner/idle ancestor gas/T/UEINT/source/fraction payloads, exact
native/replayed gas/T/UEINT/reaction/heat/face/budget outputs and independent
composite gas/species/source/numerical inventory at the unchanged `3e-14` gate (maximum residual `4.081e-16`).
Late remote-NaN, partial allocation, idle-fraction and clock mutations reject
before physics. A native-admitted provider refuses exactly once in the first
child after completed root chemistry; serial and both MPI layouts preserve all
accepted caller fields, fluxes and budgets without retry.

These are focused library diagnostics. Molecular transport is disabled in this
root fixture; optional CVODE, reacting mixtures, restart/application consumers,
full configuration matrices and physical validation were not qualified here.

The first new replay test passed its serial checks but failed when comparing MPI
ancestor fraction shapes. MPI capture stored every level's vector rather than
the selected clock's ancestor prefix; capture now stores `:selected_depth`,
including the empty root prefix. Compilation attempts exposed a fixture name
collision with `MPI_SOURCE` and overlong source lines; these were corrected.
A preflight review also separated the serial per-record reconstruction flag
from the final window acceptance flag, so a valid early record cannot conceal
late rejection. Earlier failing attempts remain in the artifact directory and
are not counted as passed checks.

Three preceding diagnostics were freshly linked and passed before these root
replay edits: zero/material ancestor capture and the unchanged external-stage
cache regression. Their exact earlier source/binary identities are retained in
`preceding-ancestor-qualification.json`; they are separate from final root replay
qualification.


Final source/library/test hashes are in `linked-identity-v5.json`; results are in
`qualification.json`, `native-replay-v5.xml`, `controls-v5.xml` and
`external-stage-final-result.json`. Full verbose replay output is retained in
`native-replay-final.log` with its separate JUnit result; this repeated execution
is one test case, not an additional qualification count.

Reproduction from the repository root:

```sh
cmake --build build/native-mol-reflux-system-clean --target \
  test_mpi_eb_native_mol_replay_3d test_mpi_eb_tree_full_3d \
  test_mpi_eb_tree_hydro_3d test_mpi_eb_held_reset_ledger_3d \
  pelef_selected_reactive_eb_3d_runtime -j 4
ctest --test-dir build/native-mol-reflux-system-clean -V \
  -R '^unit_mpi_eb_native_mol_replay_3d_3$'
env -u LD_LIBRARY_PATH -u LIBRARY_PATH -u DISPLAY \
  HWLOC_COMPONENTS=-gl,-x86 OMPI_MCA_hwloc_base_binding_policy=none \
  ctest --test-dir build/native-mol-reflux-system-clean --output-on-failure \
  -R '^(unit_mpi_eb_tree_mol_3d_3|unit_mpi_eb_tree_mol_reflux_3d_3|unit_mpi_eb_held_reset_ledger_3d_3)$'
python3 -u build/native-mol-parent-replay-20261011/run-external-regression.py
```
