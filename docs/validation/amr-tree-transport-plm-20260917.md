# Limited tree coarse donors for transport — 2026-09-17

## Source and behavior

Base: `653321e15b8a35440ac80de2ea340648bfd92957`, local `main`.
The transport API gains optional `coarse_limiter` and passes it with species
count and viscous/moving-wall data to the existing limited tree ghost helper.
The combined reaction/transport/hydro step now selects these donors for both
transport half-steps under characteristic-PLM reconstruction. Absent limiter or
PCM reconstruction preserves the previous donor path.

Verified source SHA-256:

- Transport: `1d9b0e3395da1013ff91cfea819553b9f8c46bd61c451ba99a1937dfba7bdd6f`
- Split: `92cfb0c9c5426f059d3034493df789ffb19b9bd4e97124d5543186592b5cd17a`
- Transport unit: `50333c22419c926c0833ca24089698d3a4ced97217a99749e4123c9699aa1e51`
- Split unit: `0507283f56d8445a970eca5b4826a57886f7f6f920de91c7892d33a37397a8f6`

The first two files are under `src/amr/amr_patch_tree_*_3d_mod.F90`; the unit
files are under `tests/unit/test_amr_patch_tree_*_3d.F90`.

## Verification commands

GNU Fortran 13.3.0 Debug with bounds checks/backtraces/FPE traps; MPI disabled.
Existing build: `build/selected-amr-tree-debug`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/selected-amr-tree-debug --target test_amr_patch_tree_transport_3d test_amr_patch_tree_split_3d pelef_selected_amr_tree_reactive_3d_runtime -j 6
ctest --test-dir build/selected-amr-tree-debug -R '^unit_amr_patch_tree_(transport|split)_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Fixed/selected runtime builds succeeded. Focused CTest: **2/2 passed**,
65.60 s total (transport 20.33 s, split 45.25 s), no failures or skips. Provenance,
project contract and whitespace checks passed. No Python source changed; Python
syntax verification was not applicable.

The transport unit runs both PCM and limited donors for periodic, outflow,
moving no-slip and isothermal walls on root-only and four-level trees. It checks
root serial equality, composite conservation, species closure, positive recovered
temperatures, timestep bounds and rollback after a rejected step. Nontrivial
coarse/fine transport differs from PCM in each boundary case. The split unit
exercises characteristic-PLM reaction/transport/hydro with explicit and implicit
chemistry, Euler/element boundary budgets, synchronization and rollback; its
transport-disabled PCM equality checks remain.

Logs: `build/amr-tree-transport-plm-checks/`.
Executables: `build/selected-amr-tree-debug/tests/test_amr_patch_tree_{transport,split}_3d`.

## Limits

This is focused serial implementation evidence. It does not establish global
second-order AMR convergence, an MPI/compiler/configuration matrix, external PeleC
parity, physical validation or scaling. The generic component limiter is not a
nonlinear EOS admissibility proof; reactive transport still recovers and checks
ghost primitives. Subcycle scheduling and tree MPI remain open.
