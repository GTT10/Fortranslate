# Multilevel EB flux synchronization — 2026-09-18

Implementation: `e53c0dd027dbc91f84ebd90744a2ca8a78f6a514`, based on
`9d0cdaf`. Adds `synchronize_reactive_eb_tree_flux_3d` to the core EB tree.

## Contract

Inputs are valid states already advanced with each node's raw time-integrated
open-area fluxes. The routine reconciles same-level raw faces by integer-origin
ownership, including peers with different parents. It then traverses deepest
levels first, corrects exposed cells by effective-minus-raw flux divergence,
and restricts effective fluxes by open area to the parent. Child boundary
transfers also replace the corresponding face of any neighboring parent-level
patch. This preserves conservation when only one side of a parent-level
interface is refined further. Corrections include patch and physical boundary
faces; child-covered cells are masked. Bottom-up state restriction follows.

The tree and effective flux array publish together only on success. Shape,
finite-flux, geometry and corrected-state EOS failures preserve both inputs.
This routine does not advance hydro, transport, chemistry or wall sources;
callers must account for physical exchange using the returned effective root
flux. It is a serial synchronization primitive, not a multilevel application,
subcycling algorithm, MPI qualification or physical validation.

The peer reconciliation API can now omit both state and temperature for a
flux-only transaction. Existing state-correcting calls retain their behavior.
The existing reflux API retains its default interior-interface behavior; the
new optional `all_faces` mode applies complete divergence on exposed cells.

## Completed focused verification

GNU Fortran 13.3.0, `/usr/bin/gfortran`, Debug/Ninja, existing
`build/selected-tree-mpi-debug`; no compiler/MPI configuration change.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_eb_tree_flux_3d test_eb_reactive_reflux_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_eb_(tree_flux|reactive_reflux)_3d$' --output-on-failure --output-log build/eb-tree-flux-checks/ctest.log
python3 tools/check_project_contract.py --source .
git diff --check
```

Final **2/2 passed**, zero failures/skips, total 0.19 s. The new test uses
manufactured integrated transfers and independently applies their divergence
to each node before synchronization. Six fixtures combine x/y/z planar cuts
at 0.31 with four/five-node trees: two adjacent children and either one
boundary-touching grandchild or two adjacent grandchildren of different
parents. Both refinement ratios are two. Fluxes deliberately differ by node.
Composite change agrees with independently summed effective physical boundary
exchange to at most **6.5775280774929919e-15**, normalized by
`max(1,abs(initial_component))`; the asserted bound is 3e-14.
Nonfinite deep flux and finite nonphysical corrections reject with unchanged
states/temperatures and fluxes. Existing reflux and peer-state correction
checks also pass. No hydro-step or external-solver comparison was rerun.

Provenance, project contract and whitespace checks pass. No Python source
changed. Initial compilation rejected two overlong lines introduced by the
public validator rename; line continuation fixed them. Logs are retained in
`build/eb-tree-flux-checks/{initial-build,build,ctest}.log`; per-fixture output
is in `build/selected-tree-mpi-debug/Testing/Temporary/LastTest.log`.
