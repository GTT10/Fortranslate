# EB PLM normal Riemann flux stage — 2026-09-25

Code `e03050a`, parent `2741246`. New
`compute_reactive_eb_plm_normal_fluxes_3d` connects the covered-cell normal PLM
predictor to the primitive PeleC Riemann entry. Inputs are three layers of
conserved/temperature/volume-fraction halo, target geometry, interval and
cell flattening. Per-direction outputs retain left/right primitive states,
left/right traced internal-energy density, intensive flux and Godunov face
pressure/normal velocity for subsequent transverse corrections.

Only face-normal pencils undergo EOS recovery; covered cells and unused
edge/corner state slots are skipped. Open faces require two fluid donors.
Closed faces have zero outputs. Directional rotation preserves physical
momentum components. Outputs publish together only after every face succeeds.
The stage does not multiply by face aperture or area and does not yet apply
transverse corrections, artificial viscosity or centroid interpolation.

This is the first Godunov stage API, not a selectable full EB time integrator.
Wide-halo application preparation, CTU/transverse integration, centroid flux
interpolation and input dispatch remain. Existing EB time stepping is unchanged.
Flattening is supplied by the caller, not recomputed in this stage.

## Completed checks

GNU Fortran 13.3.0, Python 3.12.3, reused Ninja Debug
`build/spray-amr-reactive-debug`. Final **1/1 CTest passed**, 1.03 s (1.08 s
total), no skips. The EB hydro regression covers x/y/z cuts with uniform and
nonuniform fields: uniform intensive fluxes agree with analytic Euler fluxes
within 3e-13, including partially open faces; Godunov pressure/velocity agree;
closed outputs are zero. Order one agrees with existing PCM PeleC fluxes within
3e-13 and differs from reconstructed nonuniform fluxes. Covered and unused
corners contain NaNs. A bad open donor rejects while retaining all prior
predictions, energies, fluxes and face data. Existing hydro checks also pass.

The first test run failed while constructing an out-of-range thermodynamic
fixture at density 1.0. It was replaced by the existing valid density scale
0.31 (nonuniform increment 0.002 per index); no solver tolerance was changed.
Both test logs are retained. Provenance, Python syntax, project contract and
whitespace checks passed; compiler warnings remain. No full matrix, application
or MPI dispatch, external full-solver parity or physical validation is claimed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_eb_hydro_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^regression_reactive_eb_hydro_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/eb-plm-normal-flux-checks/{build.log,rebuild.log,LastTest-first.log,LastTest.log}`.
