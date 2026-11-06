# Bounded eight-child family candidate

The isolated guarded-family probe now accepts optional `bounded_correction`.
Its default is false and all original curvature-screen calls retain their
behavior. No production halo option or application path is changed.

`eb_bounded_tensor_mod.F90` uses the bounded scalar half-cell differences in
three successive directions. Because nonlinear directional limiting need
not commute, it averages the six directional orders. Each order preserves
the equal-volume parent mean; averaging preserves that property and removes
a privileged axis order. This costs six directional reconstructions and is
not yet a measured solver performance result.

The candidate uses the existing volume-weighted recentering, zero-volume
placeholder handling, common family contraction, held/total-energy EOS
admission and transactional failure. It bypasses the binary curvature screen;
for this mode the legacy `smooth` output denotes no binary fallback, not a
proof that the input is physically smooth. Weighted recentering may alter
scalar bounds, so scalar boundedness is not claimed for arbitrary cut-cell
gas families.

Tests use GNU Fortran 13.3.0 with bounds/FPE checks and the existing fixed-core
archive. The tensor checker covers 20 seeded random scalar fields under all
48 signed-axis maps: 960 families. Maximum mean and covariance residuals are
both 2.220446049250313e-16. This tests coordinate permutations/reflections,
not arbitrary-angle isotropy or Euler/chemistry evolution.

The guarded probe retains its previous 164 checks and adds gas-family EOS,
the y=3 residual-band interpolation target, nonuniform volume-weighted means,
weighted gas rotation and unchanged outputs on invalid-parent failure.
The exact terminal count and return status are authoritative in
`build/eb-bounded-family-20261023/guarded-record.json`. The scalar tensor
report is in its `tensor/report.json`; both record commands and source hashes.
The fixed source HEAD is `0a0e00ebe9138145fbb0e0e745665d0d768df009` with
pre-existing dirty work.

```sh
python3 tools/check_eb_guarded_quartic.py \
  --build-dir build/eb-eos-application-20261021/isolated-openmpi/registration \
  --work-dir build/eb-bounded-family-20261023
python3 tools/check_eb_bounded_tensor.py \
  --work-dir build/eb-bounded-family-20261023/tensor
```

Next: replay complete recorded gas stencils through this family with actual
volume weights, then run a bounded same-condition evolution only if the
replay and required EOS/conservation checks support it. MPI, production
integration, fixed-PeleC trajectory agreement and physical validation are
not established. C2-C5 remain open.

Terminal result: all 190 guarded checks passed, including the weighted gas
rotation checks. Python syntax, mechanism provenance, project contract and
`git diff --check` pass. No diagnostic job remains active.
