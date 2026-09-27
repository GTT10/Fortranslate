# Bounded four-port Multijet connection — 2026-09-27

## Identity and scope

Preparation implementation: `95299ba3`. Solver source: `c4f6defc` (Fortran
unchanged through the preparation commit). Existing clean Release build
`build/decane-simple-openmpi-release`, GNU Fortran 13.3.0, OpenMPI 4.1.6,
selected decane Simple transport, MPI/spray ON, tests/SUNDIALS OFF.
Executable SHA-256:
`31b374ef2bdc0385bb464e0313a50baaa1548cd0a6692d01c7d3abc0c63bcdda`.
Pinned PeleC: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
The four source files in `references/spray_multijet_case.json` were compared
byte-for-byte with `git show` at that commit before recording their hashes.

The [case guide](../../cases/spray_multijet_3d/README.md) describes the coarse
mesh, deterministic sampling and property restrictions. In particular, the
source NC10H22 input carries a dodecane property comment; its numerical values
are preserved, not independently validated as decane data.

## Commands and completed checks

```sh
python3 tools/prepare_multijet_reference.py --bundle build/decane-3sp-simple-reference/decane_3sp.json --reference-source build/pelec-pinned-hit-source --output build/decane-multijet-connection2
# cwd: build/decane-multijet-connection2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none /usr/bin/mpiexec -n 1 ../decane-simple-openmpi-release/pelef_mpi_reactive_eb_tree_3d_selected_spray tree.nml spray.nml > fortran.log 2>&1
# cwd: repository root
python3 tools/check_jet_inventory.py --run-dir build/decane-multijet-connection2 --output build/decane-multijet-connection2/inventory-check.json
python3 -m py_compile tools/prepare_multijet_reference.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The preparer ran the mechanism provenance preflight and regenerated the pinned
Simple coefficients for comparison. All commands above completed successfully.
One bounded application run passed: two SDC steps, eight injection events,
final time `2e-7 s`, 128 gas cells and 32 retained parcels, no film.

| Inventory | Result |
| --- | ---: |
| Expected injected fuel | 1.8399999999999998e-9 kg |
| Gas fuel | 9.568729833264747e-14 kg |
| Liquid fuel | 1.8399043127016672e-9 kg |
| Reconstructed fuel relative error | 0 at reported precision |
| Budget injected-mass relative error | 1.123886701838353e-16 |
| Energy budget residual | 5.496471333632513e-15 J |
| Application composite inventory error | 7.4153676159871814e-16 |

`inventory-check.json` includes output hashes. A separate CSV check in
`port-check.json` found exactly material IDs 1–4, eight parcels per port and
equal retained liquid masses within `1e-22 kg`. These material IDs are assigned
from the port index by the application. Modified upstream input was rejected
before creating the output directory; log:
`build/multijet-source-rejection/rejection.log`.

The first attempt, `build/decane-multijet-connection/fortran.log`, failed at
injection-template validation: a single namelist radius value initialized only
the first array element. The preparer now writes the explicit four-value repeat.
That failed attempt is not counted as a passing run. MPI startup X authorization
warnings remain in both logs; the corrected run exits zero and reports completion.

No solver rebuild, full CTest matrix, installed-binary audit, rank-changing
restart, full-duration Multijet run or C++ trajectory comparison was performed
for this input-only addition. The separate full positive-rate single-jet run
remained active throughout this work.
