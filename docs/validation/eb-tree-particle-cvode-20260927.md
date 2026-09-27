# EB tree and particle CVODE SDC — 2026-09-27

## Implementation and identity

Implementation: `dbb023c0b23567af735401c59f81add84a1f1612`, parent
`eaf3e79674d40b9587d6904b677a50791e29bb52`. The matching implementation working
tree was built/tested before commit; this is not an independent clean-checkout
or installed-binary qualification. A CMake indentation adjustment after the
build does not change target configuration.

The shared fixed/selected serial/MPI EB tree frontend now accepts
`sdc_reactor='cvode'` in `&eb_tree` with SUNDIALS and reacting SDC. `native`
remains the default. One optional-procedure step helper forwards the selected
reactor through common/recursive gas stages and common/recursive particle
wrappers, including optional LES and film branches. The existing outer retry
and transaction boundaries remain in place. Fixed and selected targets link
their respective mechanism-compatible CVODE runtimes.

Serial and MPI tree checkpoints encode chemistry enablement in bit 0 and
CVODE selection in bit 1 of the existing second SDC control integer. Native
values remain 0/1; reacting CVODE uses 3. Both readers require the current
configuration's exact flags. Older readers reject the new value. Embedded
per-patch records remain split snapshots with native selector; the outer tree
record owns SDC reactor/history identity. MPI write-control consensus also
includes the named reactor. Particle restart validates the enclosed gas snapshot
before accepting the coupled continuation.

## Completed evidence

GNU Fortran 13.3.0, Open MPI 4.1.6, SUNDIALS 7.2.0 official static Fortran
interfaces; existing Debug `build/amr-cvode-app-debug`, MPI/spray/SUNDIALS ON,
selected full H2/O2 bundle. No compiler/backend configuration changed. All six
fixed/selected serial gas, MPI gas and MPI particle executable targets built.
Build-time MPI runtime dependency checks passed ABI 40. Existing real-comparison,
conversion and internal-procedure executable-stack warnings remain.

| Focused application test | Completed result |
| --- | --- |
| `application_eb_tree_sdc_cvode` | Fixed gas, common step: passed, 11.40 s |
| `application_eb_tree_sdc_cvode_recursive_selected` | Selected gas, recursive step: passed, 29.62 s |
| `application_eb_particle_sdc_cvode` | Fixed particles, recursive step and LES: passed, 57.91 s |

Three distinct tests passed, with no skips. Gas tests passed in the first
invocation. The first particle invocation failed only because its changed-reactor
negative check expected the gas-only diagnostic; the particle path correctly
rejected with `Gas snapshot restore failed`. The expected diagnostic was corrected,
active-chemistry/no-publication assertions were added, and only that particle
test was rerun (1/1 passed, total 57.96 s). The earlier failed attempt is retained
and is not counted as passed.

Gas checks use two-level 4³ plane geometry with covered/cut cells, molecular
transport, forcing and two time steps. Serial and MPI gas CSV bytes agree
exactly, including 2-to-4-rank and MPI-to-serial continuation. Both serial and MPI
reject changing CVODE to native. The composite inventory check uses `3e-12`.

Particle checks use two H2O parcels on fully regular nodal geometry, reacting
fixed gas, recursive SDC and Smagorinsky LES. Chemistry is active and LES changes
the gas state. Gas CSV, coupled budget CSV and ID-sorted particle rows agree
exactly after 2-to-4-rank continuation. Maximum normalized combined residual is
`2.4509245647058538e-15` against the unchanged `3e-12` criterion. Changed reactor,
changed/removed LES, and cut-cell LES cases reject. Reactor mismatch publishes
no final gas CSV.

Common particle execution, selected particle execution, enabled film branches,
other mechanisms and regridding were compiled/connected but not separately
run here. The selected-gas common and fixed-gas recursive combinations were
not separately run. No external solver parity or physical validation is claimed.

## Reproduction and artifacts

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/amr-cvode-app-debug --target pelef_reactive_eb_tree_3d pelef_mpi_reactive_eb_tree_3d pelef_reactive_eb_tree_3d_selected pelef_mpi_reactive_eb_tree_3d_selected pelef_mpi_eb_spray_3d pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/amr-cvode-app-debug -R '^application_eb_(tree|particle)_sdc_cvode' --output-on-failure
# After correcting the particle negative-check diagnostic:
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/amr-cvode-app-debug -R '^application_eb_particle_sdc_cvode$' --output-on-failure
python3 -m py_compile tests/application/check_eb_tree_sdc_cvode_app.py tests/application/check_les_particle_sdc_app.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, build, Python syntax, project contract and diff checks passed.
Mechanism source SHA:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Logs: `build/eb-tree-cvode-checks/`, including the original `initial-ctest.log`
and `initial-LastTest.log`, and passing `particle-ctest.log`/`particle-LastTest.log`.
Inputs, CSVs and checkpoints:

- `build/amr-cvode-app-debug/eb-tree-cvode/eb_tree_cvode_x381b00q/`
- `build/amr-cvode-app-debug/eb-tree-cvode-selected/eb_tree_cvode_lctqtfx6/`
- `build/amr-cvode-app-debug/eb-particle-cvode/les_3qsf8v9j/`
