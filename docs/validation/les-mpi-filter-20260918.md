# Distributed LES filtering — 2026-09-18

Implementation `aa6cfb9b912181696b69660a80ee97e1aa710389`, based on `5b9c9d5`.
`les_mpi_filter_mod::filter_mpi_les_field_3d` filters component-first fields on
periodic x slabs using the existing Cartesian ownership rule. Neighbor ranks
exchange one plane per iteration, propagating data across multiple owners as
needed. Widths larger than local or global x extent preserve periodic repeated
samples. The routine does not gather global field data. Transverse halos wrap
locally, then the shared separable filter evaluates the owner-local result.

Cell counts, component counts and stencil widths must agree collectively;
weights must agree exactly. Shape, finite/range checks and allocation acceptance
are collective before exchange. A final collective acceptance precedes output
publication, preserving caller output on rejected configuration or data.
The API accepts the normalized supplied stencils supported by the serial filter,
up to 65 weights. MPI packed counts are bounded. General Cartesian partitions,
physical boundaries and AMR remain open. Dynamic coefficient/flux stepping has
not yet been connected to this distributed filter.

## Independently completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug. Fixed tests-enabled
`build/spray-amr-reactive-debug` built `test_les_mpi_filter`; selected H2O2
tests-disabled `build/spray-amr-reactive-selected-debug` built
`pelef_selected_spray_runtime`, including the new MPI module.

Focused CTest **2/2 passed**, 0.96 seconds, no failures/skips: two and four ranks.
A 5-by-3-by-2 global field gives unequal slabs and one-cell local extents at four
ranks. Box ratios two/six, Gaussian ratio eight and identity filtering agree
exactly with the serial result. Different valid weights, different stencil
widths and a last-rank nonfinite field reject collectively with exact output
preservation. Provenance, Python syntax for the contract checker, project
contract and whitespace checks passed. No full matrix, selected runtime, large
MPI scaling, external solver parity or physical validation ran.

## Commands and artifacts

Logs: `build/les-mpi-filter-checks/{provenance,build,build-selected,ctest}.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_mpi_filter --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_mpi_filter_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
