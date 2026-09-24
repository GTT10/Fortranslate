# MPI recursive AMR LES integration — 2026-09-25

Implementation `a5b1dba`, based on `f67baa3`.
`advance_mpi_amr_tree_subcycle_les_3d` schedules ratio-squared child intervals
using the owner-local LES level API. Coarse endpoints provide temporal donors;
child integrated fluxes accumulate before interval reflux. Root physical-face
fluxes and synchronization corrections form the returned boundary budget.

State and temperature stay on their owning ranks. The driver clones inputs,
restricts and recovers temperatures, then publishes state, temperature, boundary
budget, maximum correction, maximum diffusivity and optional per-level counts
only after the complete recursive step succeeds. Collective allocation and
input checks preserve failure agreement. No whole-tree state gather is added.
Model, isotropic and Dynamic filter options forward to the level API.
The existing MPI selector already supports the corresponding scaled coarse
interval through `subcycle=.true.`; coupled application dispatch remains open.

## Completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpiexec`), Ninja Debug.
Fixed tests-enabled runtime/test and selected H2O2 tests-disabled runtime built.
Selected execution was not run. CTest `spray_les_mpi_amr_subcycle_2`:
**1/1 passed**, 63.94 seconds, total 63.98; no failures or skips.

The focused two-rank test uses three depths with refinement ratios 2 and 3,
requiring exact counts `[1,4,36]`. Smagorinsky and WALE exercise periodic PCM/MC
cases, Vreman exercises moving isothermal walls with an empty rank, and Dynamic
exercises distributed ownership, MC, moving isothermal walls and a single-axis
pressure-target outlet. These are four selected combinations, not a cross-product
boundary/model matrix. Owned state, temperature, boundary budget, maximum reflux
correction and maximum diffusivity equal serial results exactly. Nonowner fields
remain unallocated. Excessive dt and rank-disagreeing dt reject with unchanged
state, temperature, diagnostics and counts; the Dynamic case also checks
rank-disagreeing pressure targets. No tolerance changes were needed.

Mechanism provenance, Python syntax of the contract checker, project contract
and whitespace checks passed. Expected real-comparison compiler warnings occur
at collective exact-agreement checks. This is distributed recursive operator
evidence, not application execution, external solver parity, physical validation
or a full configuration matrix.

## Commands and artifacts

Logs: `build/les-mpi-amr-subcycle-checks/` (`build-fixed.log`,
`build-selected.log`, `ctest.log`). Launcher environment overrides below are
host-local workarounds, not application defaults.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_mpi_amr_subcycle --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_mpi_amr_subcycle_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
