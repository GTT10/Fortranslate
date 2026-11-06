# Bounded quartic half-cell differences

The isolated scalar candidate retains the quartic correction at the previously
unimproved fine y=3/4 band. Its measured initial y=3 velocity error is
0.0007402508 m/s versus componentwise MC 0.02330158 m/s. This is a saved-line
interpolation comparison, not an evolved transverse momentum result.

`bounded_quartic_slope` in `eb_extremum_slope_mod.F90` forms the equivalent
half-cell slope `(22*(q[1]-q[-1])-3*(q[2]-q[-2]))/32`. It caps its magnitude
with the EP slope bound and rejects a sign opposed to the centered slope.
Children are q0 minus/plus one quarter of the resulting slope. Difference
arithmetic removes the constant background before forming the correction.
This is an experimental adaptation of the EP bound, not the published
MUSCL/PPM scheme or pinned PeleC policy. It has no binary whole-family
curvature admission gate.

The existing EP implementation now shares the bound calculation. All 284
original scalar outputs are byte-identical to the saved pre-change outputs.
The new candidate passes the same analytic constant/linear/quadratic-average,
parent-mean, translated-step and tanh-bound checks. Across the full-domain
unit-amplitude sine controls, maximum child errors are 0.17908405,
0.01497511, 2.81927e-6 and 8.88223e-8 at 8/16/32/64 cells. The known coarse
8-cell regression against MC remains; passing algebraic checks does not
imply uniformly improved accuracy.

The runtime replay evaluates 72,000 lines from 64 recorded families. All
child pairs remain finite and preserve their parent means. Initial central
y-line errors at sampled fine y=1,3,5,7 are respectively 0.0008731854,
0.0007402508, 0.0004946198, 0.0001736874 m/s. The corresponding MC errors
are 0.02748608, 0.02330158, 0.01556962, 0.19091493 m/s. Only fine y=1..16
of this 32-fine-cell periodic domain was sampled, at parent z=6 and 8;
these values cannot replace full-domain convergence results.

Source identity: main HEAD `0a0e00ebe9138145fbb0e0e745665d0d768df009` with
pre-existing dirty changes. GNU Fortran 13.3.0, -O0, bounds/FPE checks.
Artifacts and hashes are in `build/eb-bounded-quartic-20261023/`, with
`original/` for EP regression, `candidate/` for scalar controls and
`runtime/` for the saved-line replay. Reproduce:

```sh
python3 tools/check_eb_extremum_slope.py \
  --work-dir build/eb-bounded-quartic-20261023/candidate --bounded-quartic
python3 tools/check_eb_extremum_runtime.py \
  --input build/eb-inflection-admission-20261023/grid16.log \
  --work-dir build/eb-bounded-quartic-20261023/runtime --bounded-quartic
```

Python syntax, mechanism provenance, project contract and whitespace checks
pass. No production halo or time evolution changed. Next: construct an
isolated eight-child family, preserve volume-weighted means and existing
EOS rollback/admission, and test axis ordering/rotation and discontinuity
controls before any short evolution comparison. C2-C5 remain open.
