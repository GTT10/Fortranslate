# Sign-safe common slope limiting — 2026-10-29

The16-cubed evaporation failure exposed gradient-product underflow in PLM
([original capture and fix](evaporation-full-reference-20260927.md)). Inspection
of the separate common AMReX-style slope module found the same predicate in
fourth-order interior and physical-boundary slopes. Its second-order boundary
branch uses the opposite product predicate: same-sign tiny gradients multiply
to zero and incorrectly erase a monotone slope.

The new regression scales a curved five-point profile by1e-200, -1e-200 and
1e200, covering orders2/4 and both boundary, adjacent and interior positions.
It compares the scaled slope divided by the scale with the ordinary profile
result. It also requires exact zero at the captured positive trace-species
minimum. The pre-change source fails the boundary scale-invariance assertion
(0 versus3.6666666666666665); the corrected source passes. The standalone
build enables invalid/divide-by-zero/overflow traps. The regression checks
homogeneity and extrema, not duplicated limiter arithmetic.

Only sign predicates change: direct comparisons replace products, retaining
zero-gradient semantics, formulas, order selection and face guards. No
species clipping, tolerance relaxation or reference substitution is included.
The change intentionally differs from finite-precision upstream product
behavior at extreme magnitudes; ordinary-range algebra is unchanged. It is
not yet shown to resolve the previously recorded C2 trajectory sensitivity.

Source snapshots, exact patch, negative/positive logs and commands are in
`build/amrex-slope-sign-20261029/`. The fresh n8 zero-flow outflow serial/MPI3 regular-AMR comparisons pass the
original fixed-reference field and inventory gates at all four times through
300ns. All solution and budget CSVs match the pre-change paired-ledger
baseline byte for byte, and serial/MPI3 outputs agree exactly. Maximum
x-momentum budget relative residual is2.25736e-12, below the unchanged1e-10
criterion. Three Debug tests pass with zero failures/skips: slope limiter,
EB-tree MOL spatial halos and EB-tree StateRedist (12.10s total).

`verify.py` independently checks the completed checker reports, source and
binary hashes, JUnit, baseline identity and serial/MPI identity; its result is
`verified.json`. The first verification script used an incorrect budget-file
suffix and failed before comparing budgets; correcting the path to
`fortran_NNN.budget.csv` changed no output or gate. This establishes bounded
regression coverage and extreme-value limiter behavior, not resolution of the
existing C2 maximum-norm sensitivity. Project contract, Python syntax and
whitespace checks pass.
