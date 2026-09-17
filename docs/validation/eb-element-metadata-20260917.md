# Mechanism element metadata and EB inventories — 2026-09-17

Implementation: `d9043e28efc67c00ea299b660bcd2b26b28d1504`, based on
`eae1ff5` on local `main`.

Complete mechanism bundles now generate public `<prefix>_nelements`,
`<prefix>_element_names` and `<prefix>_element_atoms(element,species)` constants.
Rows use sorted source element labels; columns retain source species order.
Labels and atom counts come from the pinned JSON thermo composition, without
inferring chemical formulae from species names. Names are escaped and wrapped
for Fortran; non-printable/non-ASCII element labels reject validation. New
symbols participate in the existing identifier-collision checks. CMake exports
their symbol names to configured frontends. Reaction-only legacy generation
does not require or emit this metadata. Generated reaction rates were not
hand-edited.

`reactive_eb_composition_integrals_3d` accepts that atom matrix and computes an
arbitrary-length elemental inventory from fluid volume and species amount
(`rhoY / molecular_weight`, giving kmol of atoms). It skips covered cells,
validates shapes, finite nonnegative composition/state and molecular weights,
and publishes the inventory only on success. Every species must have at least
one positive atom count. Chemical labels are not interpreted by this API.

This is the data/integration foundation for carbon-fuel EB support. Existing
H/O/N inventory routines and checkpoint schemas are unchanged; application
wiring, generic checkpoint element records and carbon-fuel runs remain open.

## Independently executed checks

GNU Fortran 13.3.0, Debug/Ninja, existing `build/selected-tree-mpi-debug`.

```sh
python3 tests/python/test_generate_mechanism_bundle.py \
  GenerateMechanismBundleTest.test_element_metadata_compiles_and_preserves_species_columns \
  GenerateMechanismBundleTest.test_element_metadata_names_are_safe_and_reserved \
  GenerateMechanismBundleTest.test_max_length_equation_emits_bounded_fortran_and_compiles
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug \
  --target test_eb_composition_integrals_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug \
  -R '^unit_eb_composition_integrals_3d$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tools/generate_elementary_mechanism.py tests/python/test_generate_mechanism_bundle.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Python **3/3 passed**, zero failures/skips, 0.528 s. A compiled executable
checks generated Ar/C/H/O row order and methane/CO2/Ar species columns.
Identifier collision and control-character labels reject. The existing
maximum-equation-length test compiles a complete generated bundle with the
new declarations and retains the line-length bound.

Fortran **1/1 passed**, zero failures/skips, 0.01 s (0.03 s total). The fixture
uses deliberately non-formula species aliases, a plane at z=0.6 in the unit
cube, and one kmol/m3 of each species. Integrated Ar/C/H/O amounts match
`[0.4,0.8,1.6,0.8]`. Covered storage contains NaNs and is ignored. Nonfinite
metadata, negative active species and missing composition reject without
changing the prior inventory.

Provenance, Python syntax, project contract and whitespace passed. Mechanism
source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Build/CTest logs: `build/eb-element-metadata-checks/{build,ctest}.log`.
No full generator suite, application regression, install audit, solver parity
or physical validation was run for this foundational change.
