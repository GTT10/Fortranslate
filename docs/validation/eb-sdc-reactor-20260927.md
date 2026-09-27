# EB SDC reactor callback propagation — 2026-09-27

## Source and behavior

Implementation: `8c2e626bf4fc92ec567cc9902f0539518de91b0c`, parent
`cb8a26b1131c856eff3c6b236684658c41241b9d`. Build/test used the matching
implementation working tree before commit, not an independent checkout audit.

Single-grid EB SDC, serial EB tree chemistry/common/recursive SDC and MPI EB
tree chemistry/common/recursive SDC now accept an optional
`forced_reactor_integrator`. The callback reaches each noncovered cell; absent
callbacks retain the native forced reactor. Chemistry-disabled stepping rejects
a supplied callback. MPI helpers and drivers check callback presence across
ranks before reaction, including ranks with no patches. Procedure identity
itself is not compared. Existing candidate publication and collective failure
handling preserve caller state and ledgers on rejection.

## Focused test

`unit_mpi_eb_tree_sdc_reactor_3d_2`: **1/1 passed**, 78.73 s, total 78.78 s;
zero failures/skips. Existing Debug build `build/spray-amr-reactive-debug`, GNU
Fortran 13.3.0, Open MPI 4.1.6, fixed mechanism, MPI/spray ON, SUNDIALS OFF.
No compiler/backend configuration change was made. Build warnings include
real comparisons/conversions and the internal-procedure executable-stack warning.

The three-level fixture has five 4³ patches and a z-normal plane with cut and
covered cells. Two ownership layouts distribute patches across two ranks or
assign every patch to rank 0 (one idle rank).

- Common stepping with the MPI native-reactor adapter matches serial native
  state, temperature and reaction history exactly. Callback count equals twice
  the number of owned noncovered cells (two SDC iterations).
- Serial recursive stepping through the adapter matches native MPI recursive
  state, temperature and reaction history exactly. Level counts are `[1,2,4]`.
- Inventory comparison uses the existing `2e-12` relative criterion. Composite
  closure residuals are `4.2614898476812615e-15` (common) and
  `3.9314133175908747e-15` (recursive), below the existing `3e-14` criterion.
- An adapter fails on rank 0's second cell and deliberately corrupts its private
  species/temperature arguments. Common and recursive rejection retain caller
  state, temperature, history and inventory diagnostics; recursion also retains
  step counts.
- Omission of the common-step callback on the last rank rejects before any
  callback executes, including the idle-rank layout.

The adapter invokes the existing native forced reactor. This is dispatch,
conservation and rollback evidence, **not EB CVODE execution or physical
validation**. The single-grid advance wrapper and recursive MPI presence
mismatch guard are compiled but not separately exercised. EB application input,
checkpoint reactor identity, particle wrappers and direct CVODE execution
remain subsequent integration work. No tolerances or golden hashes changed.

## Commands and artifacts

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_full_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_sdc_reactor_3d_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

All completed successfully. The test target was rebuilt after tightening its
callback-count assertion, before running CTest. Provenance source SHA:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Artifacts: `build/eb-sdc-reactor-checks/` (`preflight.log`, `build.log`,
`test-build.log`, `ctest.log`, `LastTest.log`, `contract.log`).
