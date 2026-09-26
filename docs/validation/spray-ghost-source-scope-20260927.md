# Ghost reaction source scope correction — 2026-09-27

Implementation `fe6dace4`, based on `bd89d0b8`, local `main`, version 0.246.0.
H2O2 provenance SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Reference and correction

The integrated ghost stage previously enabled direct coarse source routing.
Although that transfer preserved a standalone extensive inventory, it did not
establish the right accounting for a duplicate ghost representation of a
physical coarse particle. Feeding both reactions into coarse gas could double
count the exchange when integrating the populations.

The pinned PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`,
`Source/Particle.cpp::particleMKD` and `particleMK`, accumulates physical,
virtual and ghost reactions before one current-level `transferSource` call.
Pinned PelePhysics `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`,
`Source/Spray/SprayParticles.H::transferSource`, uses `SumBoundary` and adds to
the current level's source, including the requested ghost cells. Its PeleC
branch does not directly deposit this source into coarse valid cells.

The external checkout HEADs had moved, so this evidence was re-read with
`git show <pinned-hash>:<path>`, not inferred from current checkout contents.
Both pinned source files are copied into the artifact directory.

`advance_mpi_amr_ghost_halo` now disables direct coarse routing. Same-level
valid contributions are returned in `source`; unmatched contributions stay in
`remainder` for level-local halo source handling. The generic kick API retains
explicit coarse routing for other callers with appropriate accounting.
This supersedes the coarse-routing behavior described in the earlier ghost
stage/temporal-stage records. Recursive population/source scheduling is still
unimplemented; this correction is not a full coupled conservation claim.

## Completed checks

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6. Four tests passed in 2.80 s total:

- `spray_mpi_amr_step_2`: 0.86 s.
- `spray_mpi_eb_sample_2`: 0.56 s.
- `spray_mpi_amr_step_4`: 0.76 s.
- `spray_mpi_eb_sample_4`: 0.57 s.

The regular stage asserts zero coarse feedback, nonzero fine valid feedback
and a separate nonzero halo remainder. The EB temporal fixture places its
query at a pure fine ghost-cell center: all valid sources are zero and the
reaction remains in the halo. The combined valid-plus-remainder fluid-volume
inventory matches the frozen-gas reaction within the existing relative 1e-12
tolerance with denominator floor 1e-30. This is a diagnostic inventory check,
not permission to apply the remainder to coarse gas.

The initial build failed on an overlong test assertion; the line was split
and the failure log retained. Final build, provenance, Python syntax, project
contract and whitespace checks passed. No test failures or skips occurred.
No full matrix, installed audit, recursive application, external solver
numerical parity or physical validation was performed.

## Reproduction

Artifacts: `build/spray-ghost-source-scope-checks/`, including final CTest details.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_step test_spray_mpi_eb_sample --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_(amr_step|eb_sample)_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
