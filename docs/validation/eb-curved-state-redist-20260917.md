# Nonplanar 3D EB weighted StateRedist and coupled stepping — 2026-09-17

## Source and implementation

Base `e268adf9e0a61bcbdca147f1a8d7fa5b81680aea`, local `main`.
Replaces the axis-plane-only receiver selector with deterministic outward,
open-face neighbor selection. Receivers must meet the target volume fraction
on their own; regular and sufficiently large cut cells qualify. Existing overlap
weights, conservation algebra, mechanisms and tolerances are unchanged.

SHA-256:

- `src/eb/eb_reactive_redistribution_3d_mod.F90`:
  `21c5f9f65ba228e1a595d513be39286377ce82b8ed9d0969d495598f41dae3c5`
- `tests/unit/test_eb_reactive_redistribution_3d.F90`:
  `2544cb84b581f5529e59cce58fa5ff8c662e4b5f4bbec7c155149ee0f2742223`
- `tests/unit/test_reactive_eb_transport_3d.F90`:
  `a9253ad7d8bac3a2994fae74011c63b214b2633017770de0bd8419433667dc06`

## Commands and completed results

GNU Fortran 13.3.0 Debug; existing `build/selected-tree-mpi-debug` uses OpenMPI
4.1.6, but these focused tests execute serially.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_eb_reactive_redistribution_3d pelef_selected_reactive_eb_3d_runtime -j 4
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_reactive_eb_transport_3d -j 4
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_(eb_reactive_redistribution|reactive_eb_transport)_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Fixed core and selected runtime builds passed. An initial redistribution-only
check passed; final CTest **2/2 passed**, no failures/skips, 0.89 s total
(redistribution 0.01 s, transport/coupled driver 0.85 s). Provenance, project
contract and whitespace checks passed. No Python source changed.
Logs: `build/eb-curved-state-redist-checks/{build,coupled-build,ctest}.log`.

## Coverage

Three extrusion directions of a sampled circular section verify nontrivial
redistribution, weighted component conservation, convex value bounds, uniform
state preservation and zero covered output. Existing exact axis-plane weighting
and missing-receiver rejection checks still pass.

A 10x8x3 faceted-cylinder case advances the real coupled library driver at dt=1e-8
with hydro alone, hydro plus viscosity/conduction/species diffusion/barodiffusion,
and the full reaction/transport/hydro split. It checks positive density/temperature,
unchanged covered storage, mass/energy and elemental conservation at the existing
2e-9 bound, and nonzero separate transport and chemistry responses. Runtime chemistry
uses the fixed elementary mechanism; the selected runtime was built, not run here.

## Limits

One outward face-neighbor must independently meet the target. Thin channels or
clusters of small cells requiring chained/multi-cell merging still reject.
No nonplanar CLI/restart case, MPI/AMR EB, convergence, full matrix, installed audit,
external solver parity or physical validation was run. Faceted cylindrical geometry
approximates the sampled circle. This is not general curved-wall physical validation.
