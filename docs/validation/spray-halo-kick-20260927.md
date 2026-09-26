# Mixed-liquid particle kick from gas halos — 2026-09-27

Implementation `0998a9b1`, based on `cccfba9b`, local `main`, version 0.246.0.
Pinned PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Implementation

`kick_mpi_amr_mixture_halo` samples conserved gas and temperature guesses from
caller-filled halos at each particle's assigned node. After EOS recovery, the
existing mixed-liquid frozen-gas kernel performs drag, heating and evaporation.
Positions are restored, so this operation is a kick with no drift. Zero-mass
entries are retained without physical exchange. Matched halo deposition returns
the extensive reaction as source density per fluid volume. Input gas is read-only;
particle, component and source outputs publish together on collective success.
MPI ranks agree on the interval. Species/model tables and options must be
identical across ranks by caller contract.

The caller supplies the actual kick duration and stage-time halos. Drift,
breakup, child creation, dead-particle removal, population scheduling, source
transfer and the enclosing gas EOS/retry transaction are separate operations.
The API accepts the existing optional fluid fractions and EB visibility inputs;
this change does not establish complete recursive ghost evolution or application
integration.

## Completed checks

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6. `spray_mpi_amr_step_2` and `_4` passed
(2/2, 0.82 s and 0.78 s; total 1.65 s).

The new fixture places a mixed-liquid particle beyond the valid x extent of a
fine patch, on that patch's MPI owner. Uniform gas fills its one-cell halo.
A 1e-8 s kick matches an independent call to the existing frozen-gas kernel:
components and mass agree exactly, velocity within 1e-13 and temperature within
1e-10. Position remains exactly fixed. The global halo-volume integral of each
source component matches the kernel's extensive reaction within relative 1e-12
(with a 1e-30 denominator floor). An invalid assigned node rejects without
changing checked particle fields, component inventory or source arrays.
Existing composite and capped coupling assertions also execute.

Provenance, target build, Python syntax, project contract and whitespace checks
passed. No build/test failures or skips occurred. The new kick's EB visibility,
zero-mass and globally empty population paths were not independently exercised.
No full matrix, installed audit, recursive application, external solver parity
or physical validation was performed.

## Reproduction

Artifacts: `build/spray-halo-kick-checks/`, including copied CTest details.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_step --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_step_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
