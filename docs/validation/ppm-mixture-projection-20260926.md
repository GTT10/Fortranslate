# Independent species PPM projection — 2026-09-26

Implementation `b27bc45e`, parent `9e67823d`, local main, version 0.246.0.
The normal-projection oracle now includes independently varying positive species
mixtures. Although each cell sums to one, independently limited parabolas need
not produce face weights summing to one. Pinned PeleC retains these weights.
PeleF normalized the edges, reference integrals and final faces, changing species
and reference EOS properties before characteristic projection.

The optional `pelec_projection` cell/projection policy now preserves those
weights and uses the existing reconstructed-state EOS adapter without requiring
normalization. The regular 3D pressure-traced CTU route selects this policy and
uses the same adapter to form its intermediate conserved face states. Primitive
PeleC Riemann input also accepts nonunit positive weights through that adapter.
The default EOS-based reconstruction remains unchanged. Accepted cell-state
checks, species conservation gates and transactional publication are unchanged.
This is a numerical change to the PeleC pressure path, not a new acceptance
policy for cell compositions.

## Evidence

GNU Fortran/C++ 13.3.0, Python 3.12.3, existing Ninja Debug tests-enabled
`build/spray-amr-reactive-debug`. Expanded oracle on the old binary failed at
case 27 with scaled error **0.007305248303902001**. The corrected executable
passed **315 cases / 7,560 values**, maximum measured scaled error **0**, using
the unchanged `3e-12` acceptance threshold. This finite comparison does not imply
bitwise agreement for all inputs or compilers. Compiled source extraction,
NASA7 adapter and excluded internal-energy observable follow the
[normal oracle description](ppm-pinned-normal-20260926.md). Random seed remains
20260926; the added mixture consumes independent draws for six species, with
the seventh closing the cell composition.

Relevant CTests: **4/4 passed**, 15.84 s total, no failures/skips:
normal 3D (6.08 s), PPM step (3.41 s), halo (6.29 s), primitive Riemann (0.00 s).
The step fixture now varies two species independently before closing the third,
covering four-step conservation, species closure, physical/periodic transport
composition and rollback in both HLLC and PeleC paths. The halo unit retains the
signed-predictor versus accepted-cell rejection checks.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_step_3d test_reactive_ppm_normal_3d test_reactive_ppm_halo_3d test_reactive_pelec_riemann --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tools/check_pinned_ppm_normal.py --reference /mnt/sn850/Projects/gtt2/external/PeleC --executable build/spray-amr-reactive-debug/tests/test_reactive_ppm_normal_3d --work-dir build/ppm-mixture-projection-checks/oracle
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(unit_reactive_ppm_normal_3d|unit_reactive_ppm_halo_3d|unit_reactive_ppm_step_3d|unit_reactive_pelec_riemann)$' --output-on-failure
python3 -m py_compile tools/check_pinned_ppm_normal.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, Python syntax, project contract and whitespace passed. Build warnings
remain. Artifacts: `build/ppm-mixture-projection-checks/{red.log,red-oracle/,
build.log,build-final.log,comparison.log,ctest.log,unit.log,oracle/}`.
The oracle report records pinned sources, compiler command and executable hash.
PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
No selected build, MPI/application run, installed audit, full-CTU external
comparison or physical validation ran. Complete PeleC parity remains unproven.
