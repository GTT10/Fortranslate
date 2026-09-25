# Conserved source grids in PPM CTU — 2026-09-26

Implementation `df9a6c32`, parent `b58c9644`, local main, version 0.246.0.

Periodic and supplied-halo CTU flux APIs accept `conserved_source` arrays as an
alternative to paired primitive/internal-energy sources. The periodic wrapper
wraps source halos. The halo implementation converts cell-local sources using
the recovered donor primitives and passes them to final double-transverse
half-step prediction. Shapes, finiteness and pressure-policy compatibility are
checked; simultaneous conserved and primitive source inputs are rejected.
Conversion failure returns zero output fluxes.

## Verification

The oracle now compiles pinned `Hydro.H::pc_srctoprim` as well as the full CTU
assembly. Its NASA7 adapter supplies ideal-mixture pressure derivatives through
`rho*(gamma-1)` and `p/rho`. **38 cases / 109,440 final-flux values passed**,
maximum component-scaled error **2.988913928686505e-15**, bound unchanged at
`3e-12`. Ten additional cases exercise conserved sources with viscosity 0/0.1
and periodic/outflow/slip/no-slip/inflow boundaries. Source density, momentum,
total-energy and first two species rates are `(0.4,600,-800,400,-40000,20,-20)`;
other species rates are zero. Mixtures vary spatially, so converted primitive
sources vary despite constant conserved rates.

Normal and halo units passed **2/2**, 6.40 and 6.65 s, 13.10 s total. No failures
or skips. Provenance, Python syntax, project contract and whitespace passed.
GNU Fortran/C++ 13.3.0, Python 3.12.3, existing tests-enabled Ninja Debug build;
existing compiler warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_normal_3d test_reactive_ppm_halo_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tools/check_pinned_ppm_ctu.py --reference /mnt/sn850/Projects/gtt2/external/PeleC --executable build/spray-amr-reactive-debug/tests/test_reactive_ppm_normal_3d --work-dir build/ppm-conserved-source-checks/oracle --shock-viscosity --physical-boundaries --sources --conserved-sources
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(unit_reactive_ppm_normal_3d|unit_reactive_ppm_halo_3d)$' --output-on-failure
python3 -m py_compile tools/check_pinned_ppm_ctu.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/ppm-conserved-source-checks/` contains provenance, build,
comparison, CTest and unit logs. `oracle/` contains pinned sources, compiler log,
generated driver, inputs, outputs and report with source/shim/adapter/executable
hashes. Reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

This is a compiled reference-body comparison with an EOS adapter. Application
source generation, time integration and MPI source exchange remain open. It is
not a full external executable or physical validation.
