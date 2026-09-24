# Two-level AMR PPM hydro — 2026-09-25

Source: base `43ea448`, implementation `62f7667`. GNU Fortran 13.3.0,
Open MPI 4.1.6, Ninja Debug, Python 3.12.3. Reused fixed tests-enabled
`build/spray-amr-reactive-debug` and selected H2/O2 tests-disabled
`build/spray-amr-reactive-selected-debug` separately.
Artifacts: `build/ppm-two-level-hydro-checks/`.

The serial two-level hydro library accepts `characteristic_ppm`. Coarse and
fine steps use the existing conservative CTU operator; fine substeps fill four
halo layers with time-interpolated, limited linear coarse donors. Physical
coarse padding also expands to four layers. Periodic donor indices wrap across
the coarse domain, and the donor validator accepts the four-layer stencil.
Existing interface flux accumulation, reflux, restriction, EOS recovery and
transactional commit remain in use. The existing reaction split forwards PPM.
Application input, regridding and checkpoint acceptance are a separate follow-up.

## Checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target test_amr_reactive_chemistry_3d --parallel 2`
  and `cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_amr_reactive_3d_runtime --parallel 2`: passed.
- Initial compilation rejected a long source line and incorrect optional flux
  argument names. Both were corrected; the initial log is retained.
- Two initial test attempts rejected PPM ghost interpolation. A temporary
  CFL-based fixture timestep did not resolve the failure. Inspection identified
  the two-layer donor bounds and unwrapped periodic coarse indices; these were
  corrected and the original fixture timestep restored. Failure logs are retained.
- A third attempt passed periodic conservation but failed the full-refinement
  exact-state oracle. The fine CTU update now matches the regular-grid
  operation order and uses the same input temperature recovery path, avoiding
  a redundant pre-recovery. The exact comparison was preserved.
- A fourth attempt passed the full-refinement state oracle but rejected its
  two-pressure-outlet fixture. A focused diagnostic reproduced rejection in
  the regular-grid reference: composing four-layer pressure extrapolation at
  the corner is inadmissible for this 5% density perturbation. The fixture now
  requires exact rollback in both paths, then checks successful full-refinement
  pressure parity with a single prescribed-pressure face. Other modes retain
  both prescribed faces; no tolerance or pressure-ghost algorithm changed.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset,
  `HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`:
  `ctest --test-dir build/spray-amr-reactive-debug -R '^unit_amr_(reactive_chemistry_3d|two_level_ppm_periodic_3d)$' --output-on-failure`:
  **2/2 passed**, 115.86 s and 1.25 s, 117.16 s total; no final failures/skips.
  Checks cover nonuniform periodic conservation and CFL rejection; six touched
  boundary orientations, complete refinement and a ratio-three corner;
  physical/pressure surface conservation; exact complete-refinement parity;
  reaction/hydro composition, elemental budgets and failure rollback.
  PCM/PLM variants in the existing chemistry test also passed.
- `python3 -m py_compile tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, `git diff --check`: passed.

Selected runtime execution, MPI two-level PPM, application restart, convergence,
external reference parity and physical validation were not run. Compiler
warnings remain. Reference primitive-pressure tracing fidelity remains open.
