# MOL Lax–Friedrichs application connection — 2026-09-27

## Change and boundary

Implementation: `6aee94fbf3b914711b004626c52624a9d39238fe`, based on
`485c4d8bd7d07ccfc3b21e88fad2c6f77f08e592`. Builds/checks were run with
this implementation uncommitted. The application accepts
`mol_mode=.true.,mol_lax_friedrichs=.true.` for MPI SDC. The optional flag
passes through common/recursive gas and particle wrappers, MPI spatial rates,
and EB face evaluation to the previously checked LF kernel. MPI controls
include the flux choice. The default remains the acoustic mixture solver.
The optional `PELEF_EB_TREE_MOL_LF_V1` envelope rejects a changed flux choice
on restart; default checkpoints retain their existing envelope.

This connection is not full reference parity. Nested particle LF positivity
is unresolved, as recorded below. Film/reacting LF cases, serial application
MOL, rotating frames and extra passive/linear fields were not qualified here.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, existing Debug build
`build/heptane-simple-mpi-debug`, selected heptane Simple bundle, MPI/spray on,
SUNDIALS/tests off. Build commands (each with `env -u LD_LIBRARY_PATH
-u LIBRARY_PATH`):

```sh
cmake --build build/heptane-simple-mpi-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
cmake --build build/heptane-simple-mpi-debug --target pelef_mpi_reactive_eb_tree_3d_selected --parallel 2
python3 build/run_mol_lf_gas.py
python3 build/run_mol_lf_single.py
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
python3 -m py_compile build/run_mol_lf_app.py build/run_mol_lf_small.py build/run_mol_lf_gas.py build/run_mol_lf_single.py
git diff --check
```

Both builds and all listed final checks pass. Build logs:
`build/mol-lf-app-build.log`, `build/mol-lf-gas-build.log`.
Application drivers unset DISPLAY and both library-path variables, and set
`HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`.

- Gas: seven launches, exact final CSV on 1→2 rank restart, common/recursive
  single-level equality, bidirectional LF/acoustic restart mismatch rejection.
  `build/mol-lf-gas/result.json`; binary SHA256
  `3ccd7a3f2076f46527e5a757592e67f346b519c674383f0725dd1f464fc887df`.
- Single-level particles: seven launches, two 100 ns steps, four parcels,
  exact gas CSV and sorted particle rows on 1→2 rank restart, common and
  recursive stepping, bidirectional flux mismatch rejection. Maximum absolute
  combined ledger residual `8.131516293641283e-20`, below unchanged `1e-12`
  fixture gate. `build/mol-lf-single/result.json`; binary SHA256
  `a7d52115c858e1422fbeff8f82269c02b504befd15e70bdf48d2e65dfa1461d4`.

## Failed nested particle checks

`python3 build/run_mol_lf_app.py` fails the first 100 ns root step of the
previous two-level heated-wall particle fixture. Read-only GDB diagnosis
(`build/mol-lf-debug.gdb`, `build/mol-lf-app/lf-debug.log`) identifies the
second 50 ns fine substep: node 2, cell (2,4,4), initially zero fuel, receives
fuel rate `-0.025486868635953324` and candidate fuel partial density
`-1.2743434317976661e-9`. Thermodynamic recovery correctly rejects it.

`python3 build/run_mol_lf_small.py` reduces the root interval to 1 ns and
final time to 2 ns. Continuous stepping completes but contains negative fuel
(minimum `-2.541414850722116e-14`); the step-one checkpoint is rejected with
`Gas snapshot write failed`. This is not a passing restart qualification.
No tolerances, checkpoint checks, or species clipping were changed.
The causal source of this nested-case negative rate is not established here.
The successful single-level checks do not resolve it.
