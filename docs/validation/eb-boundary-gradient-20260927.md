# Quadratic 3D embedded-wall gradient — 2026-09-27

## Implemented boundary

`b7af009e` adds `eb_boundary_gradient_3d_mod` to the core library. It constructs
the pinned type-0 quadratic boundary stencil, including axis sorting, signed
coordinate transforms, transverse shifts and extensive area scaling. The scalar
flux API retains the reference `abs(weight)>1e-14` donor threshold and accumulation
order. It rejects active covered/nonfinite donors and publishes results only on
success. Unequal cell widths are rejected as in the pinned constructor.

This addresses a prerequisite for Spray-EB's isothermal embedded wall. The
current 3D application still assigns zero molecular flux to embedded surfaces:
the new kernel is not yet connected to transport stages, wall budgets, timestep
limits or restart. Domain-face isothermal walls are a separate existing path.
No full Spray-EB execution or physical validation is claimed.

## Evidence

Reference: PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`,
`Source/EB.cpp::pc_fill_bndry_grad_stencil_quadratic`,
`pc_apply_eb_boundry_flux_stencil` and `Source/EB.H::idxsort/argmax`.
`Source/InitEB.cpp` selects the quadratic constructor for stencil type 0;
the Spray-EB input specifies that type. The constructor normal and centroid
conventions follow `pc_fill_sv_ebg`; normalized area becomes an extensive
area-over-distance weight. The caller must gather donors using the returned
integer offset and preserve that orientation.

GNU Fortran/G++ 13.3.0, direct focused compilation:

```sh
mkdir -p build/eb-boundary-gradient-unit
/usr/bin/gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -J build/eb-boundary-gradient-unit -I build/eb-boundary-gradient-unit src/core/precision_mod.F90 src/eb/eb_boundary_gradient_3d_mod.F90 tests/unit/test_eb_boundary_gradient_3d.F90 -o build/eb-boundary-gradient-unit/check
build/eb-boundary-gradient-unit/check
python3 build/eb-boundary-gradient-unit/compare_reference.py
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile build/eb-boundary-gradient-unit/compare_reference.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Completed checks:

- One Fortran unit executable: quadratic temperature fields in 24 orientations
  reproduce the analytic normal derivative within `2e-14` absolute flux;
  anisotropic rejection preserves the previous stencil; inactive undefined
  covered donors are ignored; active covered donors preserve the previous flux.
- The comparison driver verifies `EB.cpp` and `EB.H` against `git show` at the
  pin, extracts the original constructor and sorting functions, and compiles
  them with a serial AMReX interface shim. For the same 24 orientations, all
  three offsets and 28 coefficients agree exactly at exported precision.
  Gate: `2e-14` scaled coefficient error. Report and compiler commands:
  `build/eb-boundary-gradient-unit/comparison.json`; raw outputs:
  `reference.txt`, `candidate.txt`. This is a compiled constructor comparison,
  not a complete AMReX/solver execution or independent scalar-flux comparison.
- Mechanism preflight, Python syntax, project contract and whitespace checks pass.

CTest registration is added as `unit_eb_boundary_gradient_3d`; this turn ran the
direct unit executable, not a CMake build or full CTest matrix. The existing
single-jet executable remained untouched and running.

Next integration requires gathering valid same-level/halo temperatures and
conductivity, accumulating extensive wall heat into gas and the wall ledger,
and binding wall settings to timestep/restart behavior before enabling the case.
