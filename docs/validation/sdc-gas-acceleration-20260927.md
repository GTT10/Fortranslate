# MPI EB SDC gas acceleration — 2026-09-27

Implementation `e43dd97b`, based on `4376ed09`, local `main`, version 0.246.0.
PeleC reference `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Correction

The MPI EB gas SDC common/selected-level API and spray SDC wrapper accept
optional `gas_acceleration(3)`, default zero. The required `force` argument
retains force-per-volume units. Every source evaluation computes momentum rate
`force + rho*gas_acceleration` and its work at the local stage velocity.
Acceleration must be finite, bounded by 1e6 m/s² and identical across MPI ranks.

The spray application's SDC calls now pass `external_forcing` as force density
and spray `acceleration` as both gas and particle acceleration. Previously the
application summed these unlike quantities, giving incorrect gas acceleration
at density other than one. Input/checkpoint formats are unchanged; this is a
numerical correction, and old nonzero-acceleration trajectories need not match.
The particle kick algorithm and internal-exchange ledger are unchanged.
Serial and recursive gas SDC wrappers do not expose the new optional parameter.
Recursive particle/film scheduling remains open.

## Focused verification

GNU Fortran 13.3.0, Ninja Debug, OpenMPI/OpenRTE 4.1.6, fixed-mechanism
`build/spray-amr-reactive-debug`, `/usr/bin/mpifort`, `/usr/bin/mpiexec`.

The existing force test additionally uses particle-free uniform N2 at 1000 K,
densities 2 and 4 kg/m³, velocity [1,0,3] m/s, acceleration [1000,0,3000] m/s²,
force density [20,0,40] N/m³, dt=1e-8 s and two SDC iterations. Outflow faces
and velocity tangent to the planar EB isolate uniform acceleration.
It compares total external impulse/work to M*dv and M*(v+dv/2)·dv within
3e-13 relative, cell velocity within 1e-12 m/s, temperature within 1e-8 K,
and verifies zero particle exchange. Existing ballistic particle and coupled
force checks also run.

The existing application fixture uses default regular/cut gas densities
0.31/0.93 kg/m³, nonzero gas force density and phase acceleration. It compares
2/4-rank full runs, stops on two ranks and resumes on four, checks exact gas
and particle CSV restart equality and budget agreement within 3e-12 relative,
and exercises missing, overwrite, context, hydro and truncated-file rejection.

Initial compilation caught an overlong Fortran line; it was wrapped.
The first analytic fixture with periodic x/z and slip-wall y faces rejected
its step; the final fixture uses the established outflow geometry path.
This check does not qualify that mixed-boundary setup. Both failure logs are
retained. Existing compiler warnings and X authorization messages remain.
No full configuration matrix, selected-mechanism build, installed-binary audit,
external-solver parity, temporal-order study or physical validation was run.

## Completed results

Both focused CTests passed: `spray_mpi_eb_sdc_force_2` in 13.69 s and
`spray_mpi_eb_sdc_restart` in 140.49 s (2/2, 154.18 s test time).
The density-dependent analytic assertions and changed-rank restart assertions
passed. The existing coupled force fixture's normalized residual remains
2.1293427042534180e-15. Provenance, target builds, Python syntax,
project contract and whitespace checks passed.
Application artifacts:
`build/spray-amr-reactive-debug/spray-sdc-restart/eb_spray_x_zhm0ys/`.

## Reproduction

Artifacts: `build/sdc-gas-acceleration-checks/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_reactive pelef_mpi_eb_spray_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc_force_2$' --output-on-failure
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc_restart$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py tests/application/check_mpi_eb_spray.py
python3 tools/check_project_contract.py --source .
git diff --check
```
