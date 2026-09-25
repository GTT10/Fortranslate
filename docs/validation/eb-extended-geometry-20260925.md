# Extended EB geometry and outflow CTU preparation — 2026-09-25

## Source and behavior

Implementation `0fd2a51`, parent `293b059`, adds optional `halo` to
`build_configured_eb_geometry_3d`. Default zero retains the existing physical
geometry route. Positive widths extend the computational bounds at unchanged
spacing. Nodal level-set and extruded-section input values are extended by
clamping their node indices to the physical domain; their central input values
are retained exactly. Planes use the analytic builder on the extended domain,
with wholly covered/fluid physical domains continued as such. Negative widths
fail without publishing geometry.

This follows the nonperiodic domain-face extension policy of pinned AMReX
`1546c1ccbf41c704a33790f4bf445a2fcee1372e`: `Src/EB/AMReX_EB2.cpp:31-55`
(default `extend_domain_face=true`), `AMReX_EB2_Level.H:433-436`
(nonperiodic bounding box), and `AMReX_EB2_GeometryShop.H:339-354`
(clamped implicit-function sampling). It does not claim equivalent clipping
or arbitrary AMReX topology support. The reference PeleC pin remains
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

`advance_reactive_eb_outflow_ctu_3d` accepts physical state/temperature and
geometry extended by six cells. It extracts the three-cell working margin,
copies nearest physical donors into full state halos, prepares multidirectional
flattening, and calls the complete CTU update with six physical outflow faces.
Covered halo storage is zero; an open halo requiring a covered donor fails.
Redistribution, viscosity, integrated fluxes and inventory outputs are forwarded
transactionally. Optional `shock_flattening=false` supplies unity coefficients.
It is a library entry ready for application composition; the application input
selector and reaction/transport drivers are not connected yet.

## Focused verification

GNU Fortran 13.3.0 / OpenMPI 4.1.6 / Python 3.12.3, existing Ninja Debug
`build/spray-amr-reactive-debug`; no compiler/backend changes.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_eb_3d_config --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_eb_3d_config$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_eb_hydro_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^regression_reactive_eb_hydro_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 distinct tests passed**, config 0.06 s / 0.10 s total; hydro 12.60 s /
12.64 s total. Configuration tests extend three plane directions, three
extrusion directions and one oblique nodal field by six cells. Cropping back
preserves classification and volume/face fractions, cell centroids and boundary
area within 2e-12; every extended nodal value matches the clamped input exactly,
including edges and corners. Negative-width rollback also passes.

Hydro tests exercise three plane directions with all three redistribution
choices (nine nontrivial outflow updates). Volume-weighted inventory changes
match published boundary/wall changes within 2e-12 times the inventory scale.
Invalid timestep leaves output state, temperature and ledgers unchanged.
The existing prepared CTU and physical-boundary checks also run in that test.
These flow checks cover plane geometry; extended extruded/nodal geometry is
checked by the configuration test, not a new coupled flow run.

Mechanism provenance matched SHA
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Python syntax, project contract and diff checks passed. No failed/interrupted
builds or tests, no skips and no tolerance changes. Compiler warnings remain.
Artifacts: `build/eb-extended-geometry-checks/{build.log,hydro-build.log,LastTest-config.log,LastTest-hydro.log}`.
No full configuration matrix, installed-binary audit, external solver parity
or physical validation was performed. Other boundary donors and MPI/AMR
composition remain future work.
