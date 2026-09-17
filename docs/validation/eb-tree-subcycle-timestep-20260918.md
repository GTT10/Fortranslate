# Recursive EB tree automatic intervals — 2026-09-18

Implementation: `00c4db929784d90c4a1a4f3cdb2810a8dcff67db`, based on
`d41c4c7`. Optional `subcycling` on the tree timestep API converts each patch
hydro bound to the root interval by its cumulative refinement ratio; diffusion
uses the square. The application passes its stepping mode. Common stepping
retains the previous minimum. Global donor maxima and hidden parent cells
remain included; split diffusion halves conservatively use a full-interval
bound. Checked int64 domain sizes bound the real scaling factors; comparison
before multiplication avoids overflow of inactive huge timestep bounds.

## Completed focused verification

GNU Fortran 13.3.0, `/usr/bin/gfortran`, existing Debug/Ninja build
`build/selected-tree-mpi-debug` (MPI wrapper `/usr/bin/mpifort`; tests serial).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tools/check_eb_tree_application_3d.py
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_eb_tree_subcycle_transport_3d pelef_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^(unit_eb_tree_subcycle_transport_3d|regression_eb_tree_subcycle_auto_3d)$' --output-on-failure --output-log build/eb-tree-subcycle-timestep-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^regression_eb_tree_subcycle_auto_3d$' --output-on-failure --output-log build/eb-tree-subcycle-timestep-checks/ctest.log
python3 tools/check_project_contract.py --source .
git diff --check
```

Initial focused run: **2/2 passed**, zero failures/skips, total 33.31 s
(unit 6.19 s, application 27.08 s). Inspection found that combining the checker
flags generated `subcycling=.true.0` and changed the requested maximum steps
from 200 to 20. GNU accepted that logical spelling, so the run passed, but the
checker now inserts subcycling independently in the namelist header. With the
corrected input the application rerun **1/1 passed**, 27.24 s (total 27.28 s).
These are two distinct focused tests, not three independent test cases.

The unit compares independent patch bounds with explicit substep counts
1/2/4 (hydro) and 1/4/16 (diffusion), in hydro- and diffusion-limited regimes,
for three plane orientations and four/five nodes. It checks common mode,
relaxed recursive bounds and failed-estimate rollback alongside the existing
recursive transport conservation checks.

The automatic recursive application uses a first interval of
`4.170524581499190744e-10` s below its `1e-9` s cap. Its 272 exposed cells have
maximum normalized budget residual **3.18097944570051201e-15**. Restart state
and budget CSVs match exactly. Policy/mode/corrupt-record rejection preserves
output sentinels. Provenance, Python syntax, project contract and whitespace
checks pass. No full matrix, selected-mechanism, adaptive or MPI qualification
was rerun.

Logs: `build/eb-tree-subcycle-timestep-checks/{build,initial-ctest,ctest}.log`.
Artifacts: `build/selected-tree-mpi-debug/tests/eb_tree_subcycle_auto_3d/`.
