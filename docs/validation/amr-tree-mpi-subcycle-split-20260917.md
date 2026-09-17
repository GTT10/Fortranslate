# MPI recursive split composition — 2026-09-17

## Source and implementation

Base `f4c4177db7e554040e33eb024a75e909b6c39e12`, local `main`.
The MPI tree split API now accepts optional `subcycling`, default false,
with collective agreement before dispatch. It connects the existing recursive
transport and hydro APIs inside R/2 T/2 H T/2 R/2, preserving the outer candidate
transaction and exposed-cell coarse-interval chemistry. Serial algorithms,
mechanism sources, tolerances and golden hashes are unchanged.

Verified SHA-256:

- `src/parallel/mpi_amr_patch_tree_split_3d_mod.F90`:
  `d3f52a0f3f19b317aa5ad92994926fef5e5c85366515eaa0e26933b7648c6ffb`
- `tests/unit/test_mpi_amr_patch_tree_subcycle_split_3d.F90`:
  `b025ff9150d1b2fc3e9e48585605335d1f9beeb5c322fc38c4e161c43ab5f16f`

## Focused verification

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6, existing
`build/selected-tree-mpi-debug`; no compiler/backend configuration changes.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_mpi_amr_patch_tree_subcycle_split_3d test_mpi_amr_patch_tree_split_3d pelef_selected_mpi_amr_tree_reactive_3d_runtime -j 4
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_mpi_amr_patch_tree_(subcycle_split|split)_3d_[23]$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Fixed/selected runtime builds passed. CTest **4/4 passed**, 145.58 s total:
common-step 2/3 ranks 6.04/6.04 s, recursive 2/3 ranks 66.72/66.75 s.
No failures or skips. Provenance, project contract and whitespace checks passed.

The recursive test uses a nonuniform three-level H2/O2 tree with refinement
ratios 2 and 3. Modes exercise periodic PCM/explicit chemistry,
physical moving-isothermal-wall PLM/implicit chemistry and inflow hydro with
reaction/transport disabled. Two and three ranks include empty owners.
The full interval is 1e-10 seconds. Every owned conserved state and temperature,
boundary-change vector, maximum reflux correction and minimum transport limiter
must match serial exactly. A bad hydro solver after preceding reaction/transport
work must restore all published fields and diagnostics. Rank-dependent operator
flags and subcycling mode must fail collectively without changes.

Existing common-step split tests cover default-mode compatibility. The selected
runtime library is built; the executed numerical fixture uses fixed elementary
H2/O2. Logs are in `build/amr-tree-mpi-subcycle-split-checks/{build,ctest}.log`.
No Python source changes. Existing compiler warnings remain.

## Limits

Application namelist dispatch, scaled timestep selection and mode-aware checkpoint
identity remain to be connected. This is focused implementation parity, not a
full matrix, installed-binary audit, production scaling, external-solver parity
or physical validation. Coefficient identity remains the caller's contract.
