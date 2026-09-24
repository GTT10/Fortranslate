# EB tree Mach inlet: halo and common hydro — 2026-09-25

Implementation `a852ef6`, base `9bd7ac4`.

`fill_reactive_eb_tree_halo_3d` and `advance_reactive_eb_tree_hydro_3d` now accept
optional six-face `inflow_mach`, `inflow_temperature` and `inflow_pressure`.
Negative Mach retains extrapolated outflow. Mach 0..100 activates an inward
normal inlet using interior composition and EOS reconstruction at the supplied
positive temperature. Zero pressure extrapolates interior pressure; positive
pressure prescribes it (maximum 1e15 Pa). Internal peer/coarse donors retain
their original sampling. Closed EB apertures do not evaluate inlet states.
Finite/range controls are checked even when a face is covered. Candidate halos,
tree state and boundary/wall ledgers publish only after successful completion.

This connects the previously implemented pinned nozzle state kernel. Reference:
PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`,
`build/nscbc-reference/nozzle-prob.H` SHA256 reverified as
`8e3d840a7857ec4d765019aad00369ea6fbe9aef5731fc3970c79c29e97c0bcb`.
See [state-kernel provenance](mach-inlet-state-20260925.md).

## Completed focused check

GNU Fortran 13.3.0, existing Ninja Debug fixed-mechanism build. Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_hydro_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_hydro_3d$' --output-on-failure --output-log build/eb-mach-hydro-checks/ctest.log
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 passed**, test 1.11 s, total 1.15 s. The existing hydro/conservation test
now also covers Mach 0.2, 400 K in all Cartesian orientations, extrapolated and
120000 Pa prescribed pressure, x/y/z plane cuts, and four/five-patch trees with
two refinement depths. Root ghosts match the shared state kernel exactly.
Hydro advances 1e-8 s and retains the existing composite mass/momentum/energy/
species budget tolerance 3e-14, including exposed-wall and reconciled physical
boundary fluxes. NaN temperature and missing active-inlet temperature reject
without publishing state, halos or diagnostics. Original outflow and invalid
solver/deep-state checks remain. No tolerances changed.

Build, provenance, Python syntax, project contract and whitespace checks passed;
no failures or skips. Compiler warnings remain. Artifacts:
`build/eb-mach-hydro-checks/{build.log,ctest.log,LastTest.log}`.

This is a serial library common-Euler-step increment. Recursive subcycling,
molecular transport, full split, MPI, inlet-aware timestep selection and CLI/
checkpoint input wiring are not yet connected for this EB tree mode. No selected
mechanism execution, full matrix, external solver run or nozzle validation was
performed. Existing application behavior is unchanged because the new arguments
default to outflow.
