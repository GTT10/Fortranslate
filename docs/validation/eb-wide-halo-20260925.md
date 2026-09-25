# Wide serial/MPI EB tree halos — 2026-09-25

Code `cd92508`, parent `df8a6bf`, extends serial and MPI EB tree halo entry
points with optional `halo_width` and target-resolution volume fractions.
Default width one retains existing face/full-halo behavior. Widths above one
require full-halo output and explicit target metrics, including edges/corners;
no coarse-cell metric approximation is accepted for the wide path. Central
fractions must match the target patch exactly. Physical donors are clamped to
the domain, peers take priority over ancestors, and ancestor state/temperature
interpolation follows the existing subcycling controls.

Serial buffers have shape N+2*width, representing indices 1-width:N+width.
MPI buffers use those actual spatial bounds and remain owner-local. MPI ranks
agree on width and metric presence before exchange, collectively validate local
metric fields, and publish outputs only after success. A requested open wide
halo with a covered state donor fails. Geometry/state inputs are not modified.
The underlying reference pin remains PeleC
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`; this is exchange infrastructure,
not new external-solver numerical parity evidence.

## Verification

GNU Fortran 13.3.0 / OpenMPI 4.1.6 / Python 3.12.3, existing Ninja Debug
`build/spray-amr-reactive-debug`; no compiler/backend changes.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_halo_3d test_eb_tree_halo_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_(mpi_eb_tree_halo_3d_3|eb_tree_halo_3d)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final result: **2/2 passed**, legacy serial 0.08 s, three-rank MPI 2.00 s,
2.12 s total. The MPI fixture exercises distributed owners and all patches on
rank zero (two idle ranks). Six-layer U, temperature and volume fractions
match serial exactly for root, sibling and cross-parent patches in a three-level
tree. Explicit fine metrics are preserved in every layer. Six physical outflow
layers use the same clamped boundary donor. Ancestor fractions .25/.75 also
produce exact serial/MPI agreement. Width disagreement rejects collectively
without changing output, and a central metric mismatch rejects collectively.
Existing one-layer temporal, geometry and donor-failure checks also pass.

The first build failed on two overlong Fortran lines; those were wrapped.
The first functional checks passed (1.52 s total), then the temporal wide-halo
checks were added and the final run above passed. No tolerances or golden hashes
changed, no skips. Compiler warnings remain.
Artifacts: `build/eb-wide-halo-checks/{build.log,rebuild.log,temporal-build.log,LastTest-first.log,LastTest.log}`.
Provenance SHA matched
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`;
Python syntax, project contract and diff checks passed.

## Remaining integration

These APIs supply state halos for CTU but do not yet advance a CTU patch tree.
Applications still need target-resolution extended geometry caches and
CTU stepping/reflux dispatch. Wide Mach controls are forwarded by the API but
were not newly exercised here. No MPI application, installed-binary audit,
full configuration matrix, external solver run or physical validation was run.
