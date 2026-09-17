# Sparse MPI tree conservative remapping — 2026-09-17

## Source and scope

Base: `45aa134663eff9509c26172f4a2bed693a7f45dd`, local `main`.
Adds generic/reactive MPI tree remap APIs and fixed/selected runtime build wiring.
The serial remapper only exports its existing geometry helpers; its numerical
implementation is unchanged. Application regrid dispatch is not yet connected.

Source SHA-256:

- `src/parallel/mpi_amr_patch_tree_remap_3d_mod.F90`:
  `3ea66253cd2ec771abafafa6afa0efd56c6f78ea55dd2912d6a0c9536ebb0916`
- `src/amr/amr_patch_tree_remap_3d_mod.F90`:
  `edb326d2b1f1c5f0c479e49a187cb0236bab84d646430f3e53c052e657606dc8`
- `tests/unit/test_mpi_amr_patch_tree_remap_3d.F90`:
  `4e4d9f53e0a949f310d26c332e2468d4d94eee4ab1ae9879bcd1b93c7f993ce4`

## Completed checks

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6, existing clean system-MPI configuration
`build/selected-tree-mpi-debug`. Compiler/backend settings unchanged.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 4 --target test_mpi_amr_patch_tree_remap_3d test_amr_patch_tree_remap_3d pelef_selected_mpi_amr_tree_reactive_3d_runtime
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_(mpi_amr_patch_tree_remap_3d_[23]|amr_patch_tree_remap_3d)$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Build passed for fixed and isolated selected libraries. CTest **3/3 passed**,
0.97 s total: MPI 2 ranks 0.45 s, MPI 3 ranks 0.47 s, serial remap 0.01 s.
No failures/skips. Provenance, project contract and whitespace checks passed.
No Python sources changed; no Python syntax check required. Compiler warnings
remain. No installed-prefix audit or broad configuration matrix was run.

MPI tests cover:

- A three-level old tree and a moved/branched new tree with ratios changing
  from 2/2 to 3/2, complete coarsening, and identity geometry.
- Independent old/new owner maps, including all source payload on rank zero
  and all destination payload on the last rank, with empty ranks.
- Exact serial generic output, source immutability and owner-only allocations.
  Composite two-component inventories agree within 2e-12 absolute tolerance.
- A nonuniform H2/O2 mixture: exact serial conserved state and recovered
  temperature, including a perturbed refined cell.
- Inconsistent target owners, a NaN on one source owner and a late EOS failure:
  collective rejection preserves destination state/temperature allocations and
  values. No tolerance changes or new golden hashes were used.

Logs: `build/amr-tree-mpi-remap-checks/{build,ctest}.log`.
Binaries: `build/selected-tree-mpi-debug/tests/test_mpi_amr_patch_tree_remap_3d`
and `test_amr_patch_tree_remap_3d` in the same directory.

## Limits

Payload exchange is sparse and does not gather global fields, but old exposed-cell
geometry is replicated and intersection planning is quadratic. No scaling claim
is made. Species coefficient identity across ranks remains the caller's contract.
Application event/tag selection, topology replacement and regrid-aware restart
remain to be connected. These checks establish implementation parity and bounded
conservation, not external PeleC parity or physical validation.
