# Regular spray artificial-viscosity input and restart — 2026-09-25

## Source and scope

Implementation: `082d3fdd084e35aa2c6070aa54dd3df7a0805a25`, based on
`e984ce0`. GNU Fortran 13.3.0, Ninja Debug, Python 3.12.3. Existing fixed
`build/spray-amr-reactive-debug` has tests enabled; selected H2/O2
`build/spray-amr-reactive-selected-debug` has tests disabled.

The shared serial fixed/selected application reads `difmag`, resolves automatic
selection, forwards it through single-liquid and mixture coupled stepping, and
includes the resolved PPM coefficient in checkpoint context. The coupled API
accepts a trailing optional nonnegative coefficient and validates it before
state publication. Omission retains the existing kernel default.

Input default `-1` selects 0.1 for PPM/pelec and zero otherwise. Explicit
nonnegative input requires PPM. PPM checkpoints predating this new context are
rejected, including checkpoints whose numerical default would be equivalent.
No payload/schema change was made. Default PCM/PLM context is unchanged.

## Completed checks

Commands from repository root:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_difmag_restart|spray_ppm_mach_restart)$' --output-on-failure
python3 tools/check_spray_boundaries.py --ppm --inflow --mach-inlet --difmag .05 --executable build/spray-amr-reactive-selected-debug/pelef_spray_3d_selected --work-dir build/difmag-spray-input-checks/selected
python3 -m py_compile tools/check_spray_boundaries.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Fixed CTest: **2/2 passed**, new input/restart test 24.05 s, existing PPM Mach
restart test 11.95 s (36.05 s elapsed). Selected standalone checker: **PASS**.
Provenance, Python syntax, project contract and whitespace checks passed.
No failed checks in this increment. Builds retain compiler warnings.

The bounded 4-cubed Mach-inlet water-parcel fixture includes molecular transport
and Smagorinsky. At coefficient 0.05, uninterrupted and resumed final gas and
parcel CSVs are byte-identical and final history ledgers match. Boundary mass,
energy and momentum checks retain their existing tolerances. Zero coefficient
produces a different gas CSV. Automatic (-1) and explicit 0.1 resume the same
checkpoint with identical gas, parcel and history CSVs. A changed coefficient,
negative -2, NaN, and an explicit coefficient with PLM are rejected before
history output. Existing boundary/Mach restart mismatch checks also pass.

Artifacts: `build/difmag-spray-input-checks/{build.log,selected-build.log,test.log,LastTest.log,selected-test.log}`.
Fixed new fixture:
`build/spray-amr-reactive-debug/tests/spray_difmag_restart/spray-boundaries-n7kbxtfj/`.
Existing fixture:
`build/spray-amr-reactive-debug/tests/spray_ppm_mach_restart/spray-boundaries-aap3myfj/`.
Selected fixture:
`build/difmag-spray-input-checks/selected/spray-boundaries-ay3up6sc/`.

## Evidence limits and remaining work

This is serial application/restart evidence, not an external PeleC oracle or
physical spray validation. The mixture dispatch is wired and compiled but not
separately exercised here. Selected runtime uses H2/O2, not a detailed fuel
mechanism. MPI and AMR coefficient propagation and restart binding remain open.
No full matrix, installed-binary audit or MPI runtime was rerun.
