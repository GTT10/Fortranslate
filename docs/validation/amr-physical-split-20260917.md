# Serial AMR physical-boundary reaction split — 2026-09-17

Source: parent `8bf86e7` plus this change. Logs and source hashes:
`build/amr-physical-split-checks/`. Existing GNU Fortran 13.3.0 Debug /
Open MPI 4.1.6 configuration; this check exercises the serial fixed H2/O2 API.

`advance_amr_reactive_strang_3d` now forwards optional face definitions and
reservoir states through its hydro dispatch and returns the conserved boundary
change. The ledger is private until both chemistry halves and hydro succeed.
A rejected middle hydro phase restores the entire hierarchy, including changes
from the first reaction phase, and returns zero boundary and reflux diagnostics.
Chemistry-disabled calls still dispatch hydro directly. Existing periodic callers
retain the absent-option path.

The focused unit test adds reservoir-inflow/outflow cases for PCM and
characteristic PLM with Rusanov and MC on the existing 8-by-4-by-4 hierarchy.
A reservoir density 10 percent above the initial gas exercises nonzero boundary
mass flux. Checks cover Euler inventory versus boundary change, H/O/N balance
including boundary species fluxes, nonzero chemical conversion, exact agreement
with explicit reaction-half / hydro / reaction-half calls, exact chemistry-off
hydro dispatch, and rollback after an invalid inlet rejects the middle phase.

Euler balance uses relative 2e-11 against initial inventory or boundary transfer;
all momentum components share the maximum momentum scale. Elements use relative
2e-10 against initial element inventory. These are numerical conservation and
composition checks, not physical validation.

Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target test_amr_reactive_chemistry_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^unit_amr_reactive_chemistry_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Result: the focused unit test passed 1/1 in 26.12 s, with no skips. Mechanism
provenance, project contract and `git diff --check` passed.

The first CTest selector omitted the `unit_` prefix and selected no tests; it is
not counted. A new fixture initially normalized momentum error by zero initial
momentum. Its diagnostic log records an x-momentum residual of -3.97e-23 against
8.10e-8 boundary impulse, with transverse residuals below 3.81e-26. The corrected
fixture uses the boundary impulse scale and retains the relative coefficient.
No existing tolerance, numerical solver, or Python source changed.

Physical-boundary molecular transport, MPI and application input/restart remain
open. The fine patch is still strictly interior; general multipatch/multilevel
AMR and patches touching physical faces are outside this evidence.
