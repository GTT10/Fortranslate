# Periodic PCM/PLM API viscosity — 2026-09-25

Code `c590c8c`, parent `aaa87bb`. The four periodic PCM/PLM SSPRK2 entry
points, with and without returned mean fluxes, now accept optional `difmag`.
Omission resolves to the reference 0.1. Positive coefficients use the existing
periodic boundary operator, correcting both stage fluxes before divergence.
Explicit zero preserves the previous direct periodic calculation. Negative and
nonfinite coefficients reject without publishing state; flux-returning APIs
retain zero outputs on failure.

The older Strang wrapper forwards the optional coefficient and rejects invalid
values before chemistry. Full composition and the older two-level coarse-step
fallback explicitly forward their resolved zero, preventing the new default
from overriding a user's disabled correction. Raw Riemann face-flux/RHS
primitives remain unadjusted; they do not have the geometry inputs needed for
this correction.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Python 3.12.3, Ninja Debug; reused
`build/spray-amr-reactive-debug`. **1/1 CTest passed**, 7.61 s (7.65 s total),
no failures or skips. Provenance, Python syntax, project contract and whitespace
checks passed. Compiler warnings remain.

The existing PCM/PLM viscosity test retains physical/periodic boundary checks
and adds five entry points times six policies: omitted, 0.1, 0.037, zero,
negative and NaN. Positive/default states, temperatures and mean fluxes match
the boundary operator exactly. Explicit-zero states match full composition
exactly. Published mean fluxes reconstruct accepted updates within 2e-13;
invalid inputs preserve state and temperature and zero published fluxes.
Strang checks disable chemistry to isolate forwarding. The AMR zero fallback
compiled and was inspected, but its hierarchy regression was not rerun.
No full configuration matrix, external solver parity or physical validation
is claimed. Tolerances and golden hashes were not changed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_pcm_plm_viscosity_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_pcm_plm_viscosity_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/periodic-api-difmag-checks/` contains `build.log`,
`rebuild.log`, `final-build.log` and `LastTest.log`.
