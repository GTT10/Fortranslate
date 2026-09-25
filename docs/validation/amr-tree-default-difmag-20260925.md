# AMR tree reference viscosity default — 2026-09-25

Code `3c9730e`, parent `77f2cb1`. Serial common/level and MPI common/level
regular AMR tree hydro resolve omitted `difmag` to 0.1 independently of
reconstruction and solver. Recursive stepping forwards the optional argument
through these level operators. Explicit zero still disables the correction.
Collective agreement and transactional publication are unchanged.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Python 3.12.3, Ninja Debug; reused
`build/spray-amr-reactive-debug`. **2/2 two-rank CTests passed**: common-step
6.67 s, recursive 20.73 s (27.44 s total), no failures or skips. Compiler
warnings remain. Provenance, Python syntax, contract and diff checks passed.

Each fixture covers eight ownership/boundary modes with PCM/PLM and HLLC.
Odd modes omit the serial coefficient and supply explicit 0.1 to MPI; even
modes reverse that assignment. States, temperatures, boundary and reflux
diagnostics match exactly. Recursive level step counts remain [1,2,6]. Existing
failure/rollback checks run. PPM and separate reaction/transport wrappers were
not independently rerun in this increment; no full matrix or physical validation
is claimed. Older two-level, spray and EB library defaults still require review.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_hydro_3d test_mpi_amr_patch_tree_subcycle_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_amr_tree_default_(hydro|subcycle)_3d_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/amr-tree-default-difmag-checks/build.log`, `LastTest.log`.
