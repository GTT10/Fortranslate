# Selected-mechanism LES application integration — 2026-09-27

Test integration `8e73e7aa`, based on `719a0515`, local main, version 0.246.0.
The Fortran application sources are unchanged from that base.

## Scope

The existing LES application checks now accept selected-mechanism composition
input. CTest registers common/recursive gas and particle checks when the
configured bundle matches the pinned `mechanisms/h2o2_full.json` by SHA256.
Other bundles do not inherit an incompatible fixture.

The gas fixture uses H2/O2/N2 mole proportions 2/1/3.76; the particle fixture
uses nitrogen gas and two water parcels. Both use the selected ten-species
runtime, molecular transport and Smagorinsky LES, with chemistry disabled.

## Reproduction

Fresh Ninja Debug directory: `build/les-selected-integration-debug/`.
GNU Fortran 13.3.0, Open MPI 4.1.6, Python 3.12.3. Configuration enables MPI,
spray and tests, disables SUNDIALS, and sets the bundle/source to
`mechanisms/h2o2_full.json` / `mechanisms/h2o2_cantera.yaml` using absolute paths.
Compilers/launcher: `/usr/bin/gfortran`, `/usr/bin/mpifort`, `/usr/bin/mpiexec`.
Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/les-selected-integration-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_SPRAY=ON -DPELEF_ENABLE_TESTS=ON -DPELEF_ENABLE_SUNDIALS=OFF -DPELEF_MECHANISM_BUNDLE=/home/gtt13/projects/pelec/mechanisms/h2o2_full.json -DPELEF_MECHANISM_SOURCE=/home/gtt13/projects/pelec/mechanisms/h2o2_cantera.yaml
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/les-selected-integration-debug --target pelef_amr_tree_reactive_3d_selected pelef_mpi_amr_tree_reactive_3d_selected pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/les-selected-integration-debug
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/les-selected-integration-debug -R '^les_selected_(gas|particle)_app_(common|recursive)$' --output-on-failure
python3 -m py_compile tools/check_les_mpi_amr_sdc_app.py tests/application/check_les_particle_sdc_app.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/les-selected-integration-checks/` contains configure, build,
provenance and CTest logs. Per-case inputs, output CSVs, restart artifacts and
logs are under `build/les-selected-integration-debug/les-selected-*`.

## Independently completed results

All three selected executables built successfully. Mechanism provenance,
Python syntax, project contract and whitespace checks passed. CTest completed
**4/4**, no failures or skips: gas common 46.14 s, particle common 29.81 s,
gas recursive 39.60 s, particle recursive 66.87 s; total 182.47 s.
The build retains real-comparison and executable-stack warnings.

For both time-stepping modes, LES changes the gas output relative to disabled
LES. Gas serial/MPI CSVs and checkpoints interoperate exactly, including
2-to-4-rank restart. Particle gas/budget CSV bytes and sorted particle rows
match continuous execution after 2-to-4-rank restart. Changed or removed LES
settings reject restart; particle cut-cell geometry is rejected.
Normalized combined gas/particle inventory residuals are
8.4176574716329715e-16 (common) and 2.2176423928793202e-16 (recursive), below
the existing 3e-12 criterion.

`LastTest.log` is preserved in the artifact directory. Particle run directories:

- `build/les-selected-integration-debug/les-selected-particle-common/les_xnn5ovyx/`
- `build/les-selected-integration-debug/les-selected-particle-recursive/les_cghqrzu6/`

## Limits

This is selected-runtime application integration evidence. It does not cover
active chemistry, other bundles or LES models, films, regridding, a timestep
actually shortened by the LES bound, installed binaries, a full configuration
matrix, external PeleC parity or physical validation. The regular-geometry
restriction remains; cut-cell LES is not introduced.
