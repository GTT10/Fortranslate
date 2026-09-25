# EB artificial-viscosity kernel masks — 2026-09-25

## Scope and reference

Code `8ce8baf`, parent `9a0b91a`.

The existing intensive-flux correction API now accepts an optional all-or-none
group: one-cell volume-fraction halo and three open-face masks. With a positive
coefficient, nodal divergence is suppressed when the surrounding eight-cell
volume-fraction sum is below `8*(1-1e-12)`. Closed faces become zero. Covered
cell states are not read; open faces adjacent to covered cells reject without
publishing any candidate flux. Fluid/cut states and all volume fractions are
validated. Existing no-slip momentum exclusions and species-flux normalization
remain. A zero coefficient retains the existing identity policy.

Reference inspected locally at pinned PeleC
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`: `Source/Hydro.H::eb_divu` supplies
the eight-cell cutoff; `Source/Hydro.cpp::pc_adjust_fluxes_eb` applies artificial
viscosity, normalizes species on open faces, then multiplies by face area and
aperture. This API leaves area/aperture multiplication to its caller.

EB callers do not yet supply these masks. Geometry halos, EB stepping/MPI/AMR
dispatch and application policy remain outstanding. This increment is a kernel
extension, not a complete EB hydro port or external solver parity.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug, Python 3.12.3; reused
`build/spray-amr-reactive-debug`. Focused CTest **1/1 passed**, 5.63 s
(5.67 s elapsed). No failures or skips; compiler warnings remain.

The existing PPM halo/kernel fixture now checks exact regular-limit agreement,
acceptance within the near-unit cutoff, a half-volume corner suppressing one
of four face nodes (3/4 analytic correction), unchanged distant faces, a zero
covered corner state, all-covered/closed geometry, and inconsistent open-face
rollback. Existing regular viscosity, transverse and PPM checks also run.
No tolerance was loosened. This is analytic/bounded kernel evidence, not
physical validation. Mechanism provenance, Python syntax, project contract
and whitespace checks passed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_halo_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_halo_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/eb-viscosity-kernel-checks/build.log` and copied `LastTest.log`.
