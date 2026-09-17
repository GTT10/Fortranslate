# Owner-local MPI EB checkpoint saving — 2026-09-18

Implementation: `885a2c3b116169fb579aea4618f5484148a201ec`, based on
`8febfc9`. Adds `write_mpi_eb_tree_checkpoint_3d` to fixed and selected MPI
runtimes and removes the application's remaining full-state gather helper.

## Implemented behavior

The application writes only its policy/schedule/adaptation prefix to a borrowed
scratch unit. The MPI writer appends the existing tree header and validates its
agreement across ranks. Composite inventory and checkpoint-ledger consistency
are checked collectively. Private owner-local averaging verifies the existing
2e-12 parent/child state consistency condition without modifying input fields.

Each owner then writes its patch records using the existing serial EB checkpoint
serializer, including geometry, mechanism context and patch inventory checks.
Formatted records are copied through 4096-byte buffers into an owner-local
scratch stream. Only record byte counts are reduced; node-ordered file offsets
allow owners to publish their own payloads with MPI-IO. Rank zero writes the
small common header and final marker. No patch state is gathered to rank zero.
Scratch disk usage scales with locally owned serialized records; no complete
checkpoint text buffer is held in memory.

Every patch is staged and validated before opening/truncating the destination.
Preflight failures therefore preserve existing output. I/O failures after
publication starts can leave a partial destination, as with the prior staged
serial publication; this is not an atomic filesystem replacement. The V1 tree
payload and application envelopes remain readable by the serial reader.
Initialization/restart still temporarily replicates the tree; distributed
loading is separate work.

## Completed focused verification

Fresh GNU Fortran 13.3.0 Debug/Ninja build, `/usr/bin/gfortran`, Open MPI wrapper
`/usr/bin/mpifort` and launcher `/usr/bin/mpiexec`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/mpi-eb-tree-save-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=ON
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/mpi-eb-tree-save-debug --target pelef_reactive_eb_tree_3d pelef_mpi_reactive_eb_tree_3d test_selected_full_reactive_eb_tree_3d test_selected_full_mpi_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/mpi-eb-tree-save-debug -R '^(regression_mpi_eb_tree_scheduled_3d|regression_selected_mpi_eb_tree_3d)$' --output-on-failure --output-log build/mpi-eb-tree-save-checks/ctest.log
python3 -m py_compile tools/check_mpi_eb_tree_application_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, no failures/skips, total **49.28 s**:

- Fixed explicit scheduled application: **23.84 s**, tree changes 5 -> 1 -> 4
  nodes, maximum normalized inventory residual **7.05347627465518679e-15**.
- Selected implicit H2/O2 application: **25.41 s**, static five-node tree,
  maximum normalized inventory residual **3.32177427830221154e-15**.

For both fixtures, all serialized patch records through the final tree marker
match the serial save byte-for-byte. Two/three-rank final states, two-to-three-rank
continuation and both serial/MPI restart directions match serial CSV bytes.
Scheduled restart also passes after each topology event, including the root-only
stage with empty ranks. Same-rank budget bytes match; changed-rank budget and
independent fluid-volume integrals satisfy existing bounds. Corrupt-record,
policy and selected identity/composition rejection checks preserve sentinel files.
These rejection cases exercise restart validation, not injected writer I/O faults.

Provenance, configure/build, Python syntax, project contract and whitespace
checks pass. No failed build/test attempt occurred; compiler warnings remain.
No install, subcycled/adaptive save combination, I/O fault injection, scaling,
full matrix or physical validation is claimed.

Logs: `build/mpi-eb-tree-save-checks/{configure,build,ctest}.log`.
Runtime artifacts: `build/mpi-eb-tree-save-debug/tests/{mpi_eb_tree_scheduled_3d,selected_mpi_eb_tree_3d}/`.
