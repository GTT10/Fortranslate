# Regular-geometry MPI particle LES SDC — 2026-09-27

Implementation `53e8e39d`, based on `d91427c1`, local main, version 0.246.0.

## Change

`compute_mpi_regular_eb_les_sdc_3d` adapts the existing MPI regular LES RHS to
the EB SDC diffusion interface. It collectively validates geometry and rejects
any non-regular cell or partial face. It packs separate EB parent state and
temperature endpoints into the LES temporal-donor format, selecting the active
ancestor fractions. Instantaneous face densities become interval-integrated
fluxes with the EB convention of zero-based normal indices. Outputs are
published only after collective success; intermediate allocations are checked
collectively. The caller supplies a stable timestep and matching root spacing.

Common and recursive MPI particle SDC now forward an optional diffusion callback
to the gas driver. Existing particle source iteration, synchronization and atomic
publication remain in place. This adds regular-geometry LES coupling, not a
cut-cell LES model. Application input routing is still separate work.

## Independently completed checks

Artifacts: `build/les-particle-sdc-checks/` (`provenance.log`, `build.log`,
`build-test.log`, `build-final.log`, `ctest.log`, `LastTest.log`). Existing
fixed-mechanism Ninja Debug build, GNU Fortran 13.3.0, Open MPI 4.1.6
(`/usr/bin/mpiexec`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc-les(-subcycle)?_4$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, no failures/skips: common 23.44 s, recursive 55.03 s,
CTest total 78.53 s. Build, mechanism provenance, Python syntax, contract and
whitespace checks passed. Existing real-comparison/integer-division warnings
and executable-stack linker warning remain. Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

Four-rank tests use a two-level fully regular tree, two idle ranks, nonuniform
velocity, outflow faces, molecular transport and a mixed water/peroxide parcel.
Chemistry is disabled. Smagorinsky, WALE, Vreman and Dynamic each produce positive
diffusivity and a gas state different from the no-LES reference. Normalized
combined gas/particle inventory residuals are:

| Model | Common | Recursive |
| --- | ---: | ---: |
| Smagorinsky | 1.04157e-15 | 7.70670e-16 |
| WALE | 1.24605e-15 | 9.24932e-16 |
| Vreman | 6.54320e-16 | 7.11417e-16 |
| Dynamic | 1.40866e-15 | 7.66359e-16 |

All satisfy the 5e-12 normalized criterion. Replacing regular geometry with a
valid cut-cell geometry on the second diffusion evaluation rejects the step;
gas state/temperature, component masses, particle position/velocity/temperature,
last ID, boundary/wall/chemical/forcing/exchange inventories and limiter are
retained exactly.

## Limits

No particle application input/restart route, active chemistry, physical walls,
films, deeper trees, selected mechanisms, installed-binary audit, configuration
matrix, external PeleC parity or physical LES validation was run. Geometry with
cut or covered cells is deliberately rejected by this regular LES adapter.
