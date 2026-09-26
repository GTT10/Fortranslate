# Temporal EB halos for particle gas sampling — 2026-09-27

Implementation `c6a7a484`, based on `45848ce3`, local `main`, version 0.246.0.
PeleC reference `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance preflight SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Connected path

`sample_mpi_eb_tree_particle_gas` calls `fill_mpi_eb_tree_halo_3d` with complete
corner/volume-fraction halos and caller ancestor state/temperature endpoints
and fractions. It disables snapshot restriction to retain supplied level-stage
states. The existing gas halo controls boundary/inlet data and target metrics.
It then samples conserved variables and temperature guesses at assigned
particle nodes and recovers temperature from sampled conserved variables
through the EOS. Both output arrays publish together after collective success.

This connects gas-side temporal interpolation to particle-site references.
It retains the existing PCM ancestor spatial interpolation and wider-halo
metric requirements. It does not add EB centroid visibility to the trilinear
halo sampler. Arbitrary cut-cell particle evolution, source transfer and
recursive application stepping remain open. No temporal-order or upstream
trajectory-parity claim is made.

## Completed verification

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6 through `/usr/bin/mpifort` and `/usr/bin/mpiexec`.

Four CTests passed, 2.25 s final run:
`spray_mpi_eb_sample_2`, `_4`, `spray_mpi_eb_periodic_sample_2`, `_4`.
A uniform N2 fixture at 900 K uses parent endpoint conserved states U and 2*U
and unchanged fine states. A query at the first fine ghost-cell center samples
parent fractions 0, 0.25 and 1. Actual gas-generated halo values match
(1+fraction)*U within 1e-13 normalized by max(abs(U),1), and EOS temperature
remains 900 K within 1e-8 K. The querying fine owner differs from the coarse
owner; other ranks can have empty query lists. Invalid temporal fraction 1.1
and a query outside the halo retain both prior outputs. Existing sampling
and deposition checks also run.

The initial 1000 K fixture failed temperature equality: primitive conversion
returned 1000.0000000000001 K, while recovery from the 1000 K guess returned
999.99978737271476 K. This sits on the NASA polynomial branch boundary
(`temperature <= temperature_mid` selects the low branch). The isolated
interpolation fixture was moved to 900 K; no solver tolerance or coefficient
was changed. Initial failures and diagnostic output are retained. This run
does not resolve or qualify thermodynamic round trips at the branch boundary.

Provenance, target build, Python syntax, project contract and whitespace checks
passed. No build failures occurred. Existing compiler warnings remain.
No full matrix, selected-mechanism build, installed audit, recursive application
run or physical validation was performed.

## Reproduction

Artifacts: `build/spray-temporal-halo-checks/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_sample --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_(periodic_)?sample_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
