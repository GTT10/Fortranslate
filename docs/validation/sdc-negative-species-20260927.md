# Accepted SDC species state — 2026-09-27

Implementation `00d1dc339931c8f487c4cd72d965d1dc0da9077d`, based on
`66d25ff1`. The common MPI SDC nonreacting final update and redistribution
recovery, and recursive transfer recovery, now reject every negative species
partial density. Rejection is collective and leaves published state, temperature
and ledgers untouched. Intermediate predictor recovery retains its existing
roundoff handling. No species clipping, tolerance relaxation or ledger
compensation is introduced. This repairs an acceptance/checkpoint inconsistency;
it does not solve the LF centroid-interpolation positivity problem.

## Focused evidence

GNU Fortran 13.3.0, Open MPI 4.1.6; existing Debug MPI/spray build with
SUNDIALS/tests disabled. The fixed-mechanism unit is linked manually.
Builds and tests unset `LD_LIBRARY_PATH` and `LIBRARY_PATH`; MPI execution
also unsets DISPLAY and sets `HWLOC_COMPONENTS=-x86` and
`OMPI_MCA_hwloc_base_binding_policy=none`.

```sh
cmake --build build/heptane-simple-mpi-debug --target pelef_mpi_support pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
/usr/bin/mpifort -g -fcheck=all -ffree-line-length-none -I build/heptane-simple-mpi-debug/modules -I build/heptane-simple-mpi-debug/mpi_modules tests/unit/test_mpi_mol_reaction_3d.F90 build/heptane-simple-mpi-debug/libpelef_mpi_support.a build/heptane-simple-mpi-debug/libpelef_core.a -o build/mol-negative-gate/check-final
/usr/bin/mpiexec -n 2 build/mol-negative-gate/check-final
python3 build/run_mol_negative_gate_final.py
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
python3 -m py_compile build/run_mol_negative_gate_final.py
git diff --check
```

All final commands pass. The new unit adds a stage-one-only sink of `1e-8`
kg/m³/s to an initially absent species, with matching density rate. Over
`1e-10` s the trapezoidal candidate is `-5e-19` kg/m³. Both callbacks must
execute, preventing an unrelated input rejection from satisfying the test.
The pre-fix binary accepts this candidate and fails the rejection assertion
(`build/mol-negative-gate/red-valid.log`, exit 1). The final binary rejects
it, preserving all five output ledgers, gas and temperature exactly
(`build/mol-negative-gate/final.log`, exit 0). Existing one/two-iteration
homogeneous reactor comparisons and late-reactor rollback also pass.

Six selected-mechanism nested acoustic-MOL/CTU application launches pass,
including exact 1→2 rank gas/particle restart, budgets and mode mismatch
rejection. Maximum absolute budget residual `2.710505431213761e-20`, gate
unchanged at `1e-12`. Result:
`build/mol-negative-gate-acoustic-final/result.json`; binary SHA256
`228c8227f9879952dbdb6f59d62c26b9487decc57c9eaea030b82f124492454c`.
Build log: `build/mol-negative-gate/build-final.log`.

Initial test drafts omitted the required auxiliary ledger and therefore
rejected configuration rather than a negative candidate; they are not red
regressions. A first implementation also checked intermediate predictors and
failed the CTU checkpoint case. The final implementation preserves predictor
behavior and all six application launches pass. Those intermediate results
are retained in `build/mol-negative-gate/`.

## Why the upstream cleanup is not copied into the conservative update

Pinned PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8` calls
`pc_post_eb_redistribution` from `Source/Diffusion.cpp`. In `Source/EB.cpp`,
this invokes `pc_eb_clean_massfrac` on cut neighborhoods. Defaults in
`Source/Params/_cpp_parameters` are `eb_clean_massfrac=true` and threshold 0.
`Source/Utilities.cpp::clean_massfrac` clips each species density to [0,rho],
recomputes total density, and rescales momentum/energy to preserve specific
velocity and energy. This is a nonconservative state modification.

A scalar transcription applied to the diagnosed LF candidate adds
`1.2743434317976661e-9` kg/m³ of fuel, or `9.557575738482495e-20` kg in that
cell, and changes total-energy density by `-1.4988225302658975e-5` J/m³.
It therefore changes fuel and elemental inventory without physical exchange.
Inputs, outputs and exact local upstream file SHA256 values are recorded in
`build/mol-lf-cleaner-audit.json`. This is source-derived arithmetic, not an
upstream transient simulation or proof that upstream reaches the same state.

Copying this cleanup would violate this project's species/element conservation
requirement. It is not silently enabled or emulated by relaxed checkpoint
tolerances. A conservative positivity treatment for the failing nested LF case
remains open; the checked acoustic-MOL cases remain available. The port is not
claimed complete on the strength of that alternative.
