# EB CTU conservative hydro update — 2026-09-25

Code `7e4e9ba`, parent `bdd838b`. New `advance_reactive_eb_plm_ctu_3d`
connects the twelve-solve CTU stage, artificial viscosity/species-flux
normalization, face-centroid interpolation, aperture-weighted divergence and
existing embedded-wall source. It supports ordinary conservative update,
FluxRedist and StateRedist. Accepted dt-integrated intensive centroid fluxes
and boundary/wall inventories publish together with state/temperature.

The working geometry includes three surrounding cells; the returned region is
working cells 4:N-3. The normal predictor additionally needs three halo layers
outside working geometry (six beyond the returned region). Intermediate geometry
is extracted at matching resolution. Covered output storage is canonicalized
to zero state/temperature; covered donor data are not used for EOS recovery.
All computation is private until the complete update and ledgers succeed.

Omitted viscosity resolves to 0.1; explicit zero disables its correction while
retaining species-flux normalization, which the reference applies independently
of viscosity. The existing viscosity kernel performs normalization for positive
coefficients. Centroid interpolation precedes divergence; aperture and regular
face area are applied once by the existing divergence/inventory operators.

This is a prepared-halo library update. Application input, geometry/halo
preparation, flattening and physical boundary dispatch remain. In particular,
physical no-slip masks inside the working margin are not yet supported; this
entry uses no no-slip exclusions in viscosity. Chemistry/transport/spray and
MPI/AMR CTU dispatch are not connected by this increment. No complete physical
or external-solver qualification is claimed.

## Completed checks

GNU Fortran 13.3.0, Python 3.12.3, reused Ninja Debug
`build/spray-amr-reactive-debug`. **1/1 CTest passed**, 3.42 s (3.46 s total),
no failures/skips. The EB regression adds regular and x/y/z cut geometries with
compressive nonuniform flow for all three update modes. It checks nontrivial
updates, boundary-plus-wall inventory closure within 2e-12 of inventory scale,
species closure within 2e-13, exact omitted/explicit-0.1 parity and a distinct
zero-viscosity result. Failed CFL stages preserve prior state/temperature,
ledgers and integrated fluxes. Earlier normal, transverse, CTU flux and legacy
hydro checks remain. No existing tolerance or golden was changed.

Provenance, Python syntax, project contract and whitespace checks passed.
Compiler warnings remain. No multistep application, full matrix, MPI/AMR CTU,
compiled external oracle or physical validation was run.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_eb_hydro_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^regression_reactive_eb_hydro_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/eb-ctu-update-checks/{build.log,LastTest.log}`.
