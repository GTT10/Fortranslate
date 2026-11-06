# Smooth-extremum interpolation candidate — 2026-10-23

## Decision

The pinned AMReX conservative quartic arithmetic substantially reduces the
local sine-average error, but it is **not admitted to production EB halos**.
It generates negative children at a nonnegative jump. A smooth-region guard,
admissibility treatment and EB stencil availability are required before any
evolution experiment. The existing pinned/default and opt-in policies remain
unchanged; no CFD calculation was launched for this increment.

## Reference and implementation

The local `AMReX_Interp_3D_C.H` was compared byte-for-byte with the public file
at AMReX commit `1546c1ccbf41c704a33790f4bf445a2fcee1372e`, the submodule SHA
in `references/pelec_baseline.json`. `reference.json` in the artifact directory
records the URL and SHA256. The checker rejects any different header hash.

`ccquartic_interp` forms ratio-two children from a 5×5×5 coarse stencil,
applying the one-dimensional coefficients
`2 * [-3, 22, 128, -22, 3] / 256` in z, then y, then x. Each opposite child
is formed as `2*parent - low_child`, preserving the parent average.
The isolated Fortran implementation is
`tools/reference_probes/eb_quartic_stencil_probe.F90`. It is not linked into
production halo filling or exposed as an application option.

PeleC `Source/Setup.cpp` in the saved reference source selects the linear
or piecewise-constant EB interpolaters, not this quartic interpolater.
Sharing an upstream AMReX kernel therefore does **not** establish pinned
PeleC trajectory parity for a future policy using it.

## Completed checks

Artifacts: `build/eb-quartic-stencil-20261023/`.
`tools/check_eb_quartic_stencil.py` extracts the actual pinned C++ function,
provides array-index adapters, compiles it and the Fortran implementation,
and compares 119 stencils: 15 analytic polynomial cases (including mixed
3D terms), 100 seeded random stencils, three exact sine-average stencils,
and one nonnegative jump. The reference adapter exercises child parity at
one coarse origin; it does not test AMReX's general index layout or EB logic.

- Maximum C++ versus Fortran difference: **0**.
- Maximum parent-versus-eight-child mean residual: **2.221e-16**.
- Maximum error against independently integrated polynomial children:
  **1.388e-17**.

| Root cells | Limited/constant donor error | Centered linear error | Quartic error |
| --- | ---: | ---: | ---: |
| 8 | 0.7417904 | 0.0825491 | 0.0101413 |
| 16 | 0.1909149 | 0.00546732 | 0.000173687 |
| 32 | 0.0480754 | 0.000346687 | 0.00000277674 |

Values are velocity errors in m/s for the analytic fine average next to the
shear maximum. They are local interpolation estimates, not new CFD results
or formal time-evolution convergence orders. At 16 cells the centered and
quartic errors are respectively about 35 and 1099 times smaller than the
limited donor error. The centered result is an arithmetic diagnostic only;
it does not justify disabling limiting generally.

For the stencil varying as `[0,0,0,1,1]` in y and constant in x/z, quartic
children range from **-0.1484375 to +0.1484375**, while their mean remains
zero. This is unacceptable for a nonnegative species-density field despite
exact conservation. The report explicitly records
`production_candidate_admitted=false`; the negative test is a measured
candidate limitation, not a passing positivity test. Gas EOS admission,
cut-cell geometry, neighboring stencil connectivity and MPI are untested.

Reproduction:

```
python3 tools/check_eb_quartic_stencil.py \
  --header build/eb-quartic-stencil-20261023/AMReX_Interp_3D_C.H \
  --work-dir build/eb-quartic-stencil-20261023
```

`report.json` records compiler commands and source/reference/executable hashes;
`reference.txt` and `probe.txt` retain all outputs. Compilers were the current
system C++ compiler and GNU Fortran 13.3.0, with optimization disabled and
Fortran bounds/floating-point checks enabled. Mechanism provenance, Python
syntax, project contract and whitespace checks pass. This adds an isolated
candidate and rejection evidence, not an adopted numerical correction or C2
completion.
