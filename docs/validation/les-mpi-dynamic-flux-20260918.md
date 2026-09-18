# Distributed dynamic LES fluxes — 2026-09-18

Implementation `43b9c98ee4cea200e8eabcadd052cf402c8b3e6c`, based on `5893e88`.
`compute_mpi_dynamic_les_rhs_3d` connects distributed dynamic coefficients to
owner-local periodic SGS flux divergence. State, temperature guesses and
coefficients exchange two x halo layers by neighbor-plane propagation, including
one-cell local slabs. The shared serial flux kernel accepts supplied coefficient
fields and an active x interval requiring two halo cells on either side. Only
faces touching owned cells contribute to the diffusion estimate; output outside
the active interval is zero. Transverse directions remain periodic.

Each rank computes the same arithmetic face interpolation from matching halo
values; no global field is gathered. Local acceptance is reduced before RHS
publication, and the diffusion estimate is reduced by a global maximum. Invalid
local input preserves all caller RHS arrays and returns zero diffusion output.
The existing frozen-coefficient estimate is not a general nonlinear stability
proof. MPI time integration, spray/application dispatch, physical boundaries
and AMR dynamic LES remain open.

## Independently completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; unchanged configurations.
Fixed tests-enabled `build/spray-amr-reactive-debug` built `test_les_mpi_dynamic`
and `test_spray_core`. Selected H2O2 tests-disabled
`build/spray-amr-reactive-selected-debug` built `pelef_selected_spray_runtime`.

Focused CTest **4/4 passed**, 1.13 seconds, zero failures/skips: MPI dynamic
checks at two/four ranks plus existing serial dynamic flux and SSPRK2 checks.
The MPI fixture uses a 5-by-3-by-2 variable-density/temperature/velocity field,
anisotropic spacing, unequal slabs and one-cell local extents. Box and Gaussian
configurations match serial RHS and global maximum diffusion exactly. Global
component flux sums pass the unchanged `1e-12` relative gate. A last-rank invalid
EOS guess rejects with exact RHS preservation and zero diffusion output. The
previous coefficient parity/configuration rejection checks remain.
Provenance, Python syntax for the contract checker, project contract and diff
whitespace checks passed. No failed builds/tests occurred. Selected runtime was
compiled, not exercised; no full matrix, external solver parity, large-scale MPI
or physical validation ran.

## Commands and artifacts

Logs: `build/les-mpi-dynamic-flux-checks/{provenance,build,build-selected,ctest}.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_mpi_dynamic test_spray_core --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_(mpi_dynamic_[24]|dynamic_flux|dynamic_step)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
