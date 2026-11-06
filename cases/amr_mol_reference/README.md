# Spatial-gradient AMR MOL reference control

This frozen, periodic gamma-law HIT control imports the pinned PeleC initial
cell states into the production Fortran tree integrator. It uses a 4³ root,
one 4×8×8 fine patch at ratio two, and three 100 ns steps. Chemistry and
transport are disabled. It checks the regular AMR operator, not the gas
application initializer, reacting flow, EB, regridding, MPI or convergence.

Build `PeleC-HIT` from the pinned source, the existing plotfile exporter, and
the `test_amr_tree_mol_reference` target in a tests-enabled PeleF build. Then:

```sh
python3 tools/check_amr_tree_mol_reference.py \
  --pelec /absolute/path/to/PeleC-HIT \
  --reference-source /absolute/path/to/pinned/PeleC \
  --exporter /absolute/path/to/pelef_plotfile_csv \
  --probe /absolute/path/to/test_amr_tree_mol_reference \
  --work build/amr-mol-reference --case sub1_it1
```

The source must contain the baseline provenance marker. The checker verifies
the source hashes in `references/pelec_amr_mol.json`, captures native plotfiles
with a SHA manifest, and records commands, binary hashes, times and errors in
`checks.json`. To reuse native results, replace `--pelec` and
`--reference-source` with `--reference /absolute/path/to/captured/reference`.
No reference simulation is repeated in that mode.

`sub0` selects common clocks, `sub1` recursive clocks; `it1`/`it2` select one/two
MOL iterations. Omitting `--case` attempts all four and stops on the first
failing case. All four physical times within a numerically mismatched case
are retained. A nonzero exit is a failed comparison, even if both solvers ran.

Use `--keep-going` to collect all requested clock/iteration cases despite
field or inventory gate failures. Each completed case records its status;
`matrix_complete=true` means all requested cases were compared, not that
they passed. Any failed case keeps the overall status failed and exit code
nonzero. Provenance, execution and malformed-data failures still stop the
checker immediately. Neither mode changes acceptance limits.

The predeclared field tolerance is `1e-10` times the reference field maximum
(floor `1e-12`); initial agreement is eight double-precision epsilons. The
composite mass difference must be below `1e-11` relative. Actual times, grid
spacing, all cell indices, finite fields and nontrivial evolution are checked.

**Known open difference:** nonreacting two-iteration MOL uses an old-half
flux register in the pinned reference. PeleF conserves the accepted two-half
transfer separately. The recursive `sub1_it2` comparison fails these unchanged
criteria; it is not an accepted alternate baseline. See the
[dated record](../../docs/validation/amr-mol-reference-20261019.md).
