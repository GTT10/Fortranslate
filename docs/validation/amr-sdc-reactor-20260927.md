# AMR SDC reactor dispatch — 2026-09-27

## Change and identity

Implementation: `43fbbd793892e291c052efd424aeb5a3aec87d8b`, parent
`063b14845d412944d3acf09210f80418bbcbe2ec`. The focused build/test ran on the
identical implementation working tree before its commit, not an independent
clean-checkout qualification. No numerical tolerance or golden data changed.

Serial and MPI AMR chemistry helpers and common/recursive SDC drivers accept
an optional `forced_reactor_integrator`. Each selected cell receives the same
callback through level recursion. An absent callback retains native chemistry.
A callback with chemistry disabled is rejected. MPI compares callback presence
collectively, including idle ranks; callback identity itself is not compared.
Candidate state, temperature, reaction history and step diagnostics remain
transactional across local failure and the collective result.

## Focused evidence

GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/gfortran`, `/usr/bin/mpiexec`),
existing Debug build `build/spray-amr-reactive-debug`, fixed mechanism,
MPI/spray enabled and SUNDIALS disabled. No build configuration was changed.

`unit_mpi_amr_patch_tree_sdc_3d_4`: **1/1 passed**, 25.84 s
(CTest total 25.89 s), zero failures/skips. The existing three-level test
assigns three patches to ranks 0 and 1, leaving two ranks idle. Added checks:

- Serial/common and MPI/common callback execution and exact state/history and
  inventory agreement; only owners execute MPI callbacks.
- MPI recursive callback results equal serial native results exactly.
- An adapter fails on rank 0's second cell, after a successful cell, and
  deliberately overwrites its private species/temperature arguments. Common
  and recursive calls reject and retain state, temperature, history, inventory,
  diagnostics and recursive step counts.
- Missing callback on one idle rank rejects recursive stepping before any
  reactor callback executes.

The adapter invokes the existing native forced reactor. This is dispatch and
rollback evidence, **not AMR CVODE execution**, application selection, restart
policy qualification or external solver parity. Serial recursive callback
propagation is compiled; its success path is not independently exercised here.

## Commands and artifacts

Run from the repository root:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_split_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_amr_patch_tree_sdc_3d_4$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

All completed successfully. Provenance SHA:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Logs: `build/amr-sdc-reactor-checks/` (`preflight.log`, `build.log`,
`ctest.log`, `LastTest.log`, `contract.log`). Compiler warnings include real
comparisons and the internal-procedure executable-stack linker warning.
