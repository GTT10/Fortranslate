# Distributed dynamic LES coefficient fields — 2026-09-18

Implementation `669d28108448d3d446429a737b200997eae25ea8`, based on `cf81ae4`.
`les_mpi_dynamic_mod::compute_mpi_dynamic_les_coefficients_3d` computes the three
dynamic coefficients on periodic x slabs. Conserved-state filtering, EOS recovery,
centered primitive gradients, grid/test quantities, filtering, Germano contractions
and ratio-six box smoothing operate on owner-local arrays. Primitive gradients
use one-plane neighbor exchange; wider filters use the preceding distributed
filter. No global state is gathered. Nonuniform slabs and one-cell local x
extents are supported; global extents must be at least two in each direction.

Model, filter/grid ratios, species count, global cells and spacing are collectively
checked. Allocation and local numerical acceptance are reduced before subsequent
communications. Output is published by the final transactional distributed
filter only after all ranks accept. As with existing solver APIs, callers supply
the same mechanism data on all ranks; species count agreement is checked here,
not a full mechanism fingerprint. Distributed dynamic fluxes, timestep selection,
SSPRK2, spray coupling and applications remain to be connected.

## Independently completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; unchanged configurations.
Fixed tests-enabled `build/spray-amr-reactive-debug` built `test_les_mpi_dynamic`.
Selected H2O2 tests-disabled `build/spray-amr-reactive-selected-debug` built
`pelef_selected_spray_runtime`, including the new MPI module.

Focused CTest **2/2 passed**, 0.98 seconds, zero failures/skips: two/four ranks.
A 5-by-3-by-2 variable-density/temperature/velocity field with anisotropic spacing
uses unequal slabs and one-cell local extents. Box test ratio two/grid one and
Gaussian test ratio four/grid two yield exact serial coefficient equality.
Coefficients are nontrivial and inputs remain unchanged. Different grid widths,
a last-rank invalid model, and a last-rank nonfinite EOS temperature guess reject
collectively with exact coefficient preservation. Provenance, Python syntax for
the contract checker, project contract and whitespace checks passed. No full
matrix, selected runtime, external solver parity, large-scale MPI or physical
validation ran.

## Commands and artifacts

Logs: `build/les-mpi-dynamic-field-checks/{provenance,build,build-selected,ctest}.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_mpi_dynamic --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_mpi_dynamic_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
