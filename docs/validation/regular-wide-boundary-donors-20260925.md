# Wide regular-grid boundary donors — 2026-09-25

Implementation `4bc8c45`, based on `5702d67`.
`fill_reactive_boundary_ghosts_3d` now accepts optional `ghost_width` in 1..128,
with default 2. Shape checks use that width and guard dimension overflow.
Periodic and open donors support the wider coordinates. Wall mapping repeats
until the donor lies inside the domain, then applies encountered reflections
in reverse order. This preserves distinct prescribed tangential velocities at
opposing walls. The existing momentum/kinetic-energy correction, temperature
copy, inlet precedence and axis order are retained. Euler still reflects only
normal momentum; diffusive no-slip reflects all components about wall velocity.

This supplies the regular-grid donor dependency needed by Dynamic LES physical
boundaries. It does not yet connect those boundaries to the regular Dynamic
RHS, integrator or application. The AMR path already has its own wide sampler.
The helper retains its existing intent(out), zero-initialized failure contract;
it is not an atomic state-update API.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug. Fixed core/unit target and selected H2O2 spray
runtime built. Selected execution and MPI execution were not run.
`unit_amr_patch_tree_ghost_3d`: **1/1 passed**, 0.12 seconds (total 0.16).

On a four-cell-per-axis domain, the new 13-layer regular donors equal the AMR
sampler exactly for two moving walls with different speeds, slip walls, and
a wall paired with outflow. Those states also match the existing independent
analytic repeated-reflection oracle, including momentum and energy. All copied
temperatures remain 800 K. Width 0 and 129 reject. The same test retains its
default-width physical, periodic and temporal-donor comparisons. No failures,
skips or tolerance changes occurred. Mechanism provenance, Python syntax of
the contract checker, project contract and whitespace passed. Existing compiler
warnings remain in dependent modules. This is boundary-operator evidence, not
Dynamic physical-boundary integration or external solver/physical validation.

## Commands and artifacts

Logs: `build/regular-wide-ghost-checks/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_amr_patch_tree_ghost_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
ctest --test-dir build/spray-amr-reactive-debug -R '^unit_amr_patch_tree_ghost_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
