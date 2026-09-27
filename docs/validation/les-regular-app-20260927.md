# Regular gas LES SDC application integration, 2026-09-27

## Implemented boundary

Source: `dc2b26b7b994594bb6f2e79f918a91901369de26`, based on
`1f8665cf` on local `main`. Fixed and selected regular 3D gas applications
now read optional `reactive_les` input when spray/LES support is built.
Smagorinsky, WALE, Vreman and Dynamic contribute old/new SDC source rates
alongside molecular transport, forcing and native or explicitly selected
CVODE chemistry. The core application receives an optional diffusion callback;
LES model implementation remains in the optional runtime. Fixed and selected
executables link their corresponding runtime, preserving mechanism isolation.

The timestep is bounded by `0.1/(maximum_diffusivity*sum(1/h**2))`, matching
the AMR application convention. The final trapezoidal LES boundary inventory
is moved from the source ledger into the boundary ledger. The gas SDC driver
retains its existing candidate transaction. Split stepping rejects enabled
LES. Builds without spray/LES reject enabled input rather than ignoring it.

This is application integration, not physical LES validation, reference solver
parity, a full boundary/model matrix, or an installed-binary audit. Regular gas
applications do not gain checkpoint/restart here. Timestep-limit dominance was
not separately exercised by the short-duration fixture.

## Completed checks

Used existing Debug caches without changing compiler or backend configuration:
GNU Fortran 13.3.0, `/usr/bin/gfortran`, SUNDIALS 7.2.0 static Fortran modules.
`build/amr-cvode-app-debug` has MPI/spray/SUNDIALS enabled, Open MPI 4.1.6,
`/usr/bin/mpifort`, and selected `mechanisms/h2o2_full.json`.
`build/forced-cvode-debug` has MPI/spray disabled and SUNDIALS enabled.
Existing real-comparison and executable-stack warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/amr-cvode-app-debug \
  --target pelef_reactive_3d pelef_reactive_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest \
  --test-dir build/amr-cvode-app-debug -R '^application_les_regular$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/forced-cvode-debug \
  --target pelef_reactive_3d --parallel 2
python3 -m py_compile tests/application/check_les_regular_app.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, both builds, syntax, project contract and whitespace checks pass.
Pinned YAML SHA-256 remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Final CTest: **1/1 passed, 0 failed, 0 skipped**, 11.26 s (total 11.31 s).
The test uses 4³ cells, no-slip walls on x/y, periodic z, reacting full H2/O2,
molecular transport, forcing and CVODE, advancing 1e-10 s with two SDC
iterations. All four LES closures change the gas relative to disabled LES;
fixed and matched-composition selected CSVs are byte-identical. Maximum
normalized Euler residual is `2.5582913394502576e-17`, below `3e-12`.
Unknown models and split+LES input reject without publishing the final CSV.

A separate direct check ran the spray-disabled binary against the successful
Dynamic input and verified nonzero exit, `Regular LES backend is unavailable`,
and absence of `sdc_reacting.csv`. The input is under
`build/amr-cvode-app-debug/les-regular-app/les_regular_jd74vbz8/dynamic/input.nml`;
the rejection log is `build/les-regular-app-checks/unavailable/run.log`.

Development failures were corrected before the final run: an internal-procedure
interface declaration failed compilation; passing an unassociated procedure
pointer failed Debug checks (replaced by an optional-procedure dispatch helper);
the initial fixture used invalid boundary keys; and one-direction wall shear
produced zero Vreman stress. The final fixture supplies x/y wall gradients.
No acceptance tolerances or golden data changed.

Build/test logs and initial failure logs are in `build/les-regular-app-checks/`.
Final per-case inputs, CSVs and logs are under
`build/amr-cvode-app-debug/les-regular-app/les_regular_jd74vbz8/`.
