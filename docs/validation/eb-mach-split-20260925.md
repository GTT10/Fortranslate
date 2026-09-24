# EB tree Mach inlet in transport and split stepping — 2026-09-25

Implementation `42368d4`, base `5538554`. This extends the
[common hydro connection](eb-mach-hydro-20260925.md) through serial common and
recursive molecular transport, recursive hydro and the complete
reaction/transport/hydro split. Each stage/substep reconstructs inlet ghosts
from its current donor state. Active Mach inlet faces permit molecular flux;
inactive physical faces retain closed transport. Peer/ancestor donor sampling,
interval reflux, exposed-wall accounting and transactional commit remain.
The three optional six-face inputs retain their halo API semantics.

The nozzle state kernel is unchanged: pinned PeleC
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`, nozzle-prob.H SHA256
`8e3d840a7857ec4d765019aad00369ea6fbe9aef5731fc3970c79c29e97c0bcb`.
No new reference equations were introduced.

## Focused verification setup

GNU Fortran 13.3.0, Ninja Debug, fixed elementary H2/O2 mechanism.
The five-patch, two-refinement-depth plane-cut test uses low-x Mach 0.2 with
extrapolated pressure and high-y Mach 0.1 with 120000 Pa prescribed pressure;
both use 400 K. It advances 1e-8 s in common and recursive modes, with all four
chemistry/transport enable combinations and explicit/implicit chemistry.
Composite conserved budgets retain tolerance 3e-14. Explicit full splitting
must exactly match separately invoked operators in all cell states,
temperatures and diagnostics. Inlet transport energy flux must be nonzero.
Invalid inlet temperature after a private reaction stage must roll back every
state and diagnostic, alongside the existing invalid-hydro rollback check.

The first build rejected an overlong test line; it was wrapped. Initial Mach
runs passed conservation and exact cell-state comparison but rejected exact
ledger comparison: the test summed T1+H+T2 while the public split defines
H+(T1+T2). Previously closed boundaries masked this floating-point ordering
difference. The independent composition now groups its ledger terms according
to that definition; no tolerance was changed. Failure logs are retained.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_subcycle_full_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_(mach_subcycle_full|mach_common_full|subcycle_full)_3d$' --output-on-failure --output-log build/eb-mach-split-checks/ctest.log
# After correcting test ledger grouping, rerun the two failed Mach tests:
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_mach_(subcycle|common)_full_3d$' --output-on-failure --output-log build/eb-mach-split-checks/ctest.log
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final Mach rerun: **2/2 passed**, recursive **47.96 s**, common **35.43 s**,
total **83.43 s**. The inactive-inlet baseline passed in the initial run in
**48.32 s**, before the test ledger grouping correction; it was not repeated.
Build, provenance, Python syntax, project contract and whitespace passed.
No final failures or skips; compiler warnings remain. Artifacts in
`build/eb-mach-split-checks/`: build.log, build-initial-failure.log,
ctest.log, LastTest.log, ctest-initial-failure.log and LastTest-initial-failure.log.

This remains a serial library increment. MPI, inlet-aware timestep bounds and
application input/checkpoint wiring remain open. No selected-mechanism execution,
full matrix, external-solver comparison or physical nozzle validation is claimed.
