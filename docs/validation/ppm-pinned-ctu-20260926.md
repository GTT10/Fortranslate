# Compiled pinned periodic CTU assembly — 2026-09-26

Comparison implementation `6861c0b8`, parent `23cd12ac`, local main, version
0.246.0. No solver arithmetic, tolerance or golden data changed.

`tools/check_pinned_ppm_ctu.py` extracts unchanged pinned PeleC `pc_umeth_3D`
and its PPM tracing, Riemann, passive-upwind, transverse and pressure-divergence
function bodies. Thus the reference determines the twelve-solve operation order
and intermediate array indices; the checker does not rewrite the assembly in
Python or Fortran. Final `pc_norm_spec_flx` is also compiled unchanged.
`tools/reference_shims/ctu_array.hpp` supplies serial loops, bounded array views
and storage. Unsupported PLM/WENO/inflow/flattening branches terminate if entered.
The existing NASA7 SI adapter consumes the fixed Fortran thermo coefficients.
This is a compiled reference-body oracle, not an installed AMReX/PeleC executable.

## Results and scope

**6 cases / 17,280 final-flux values passed**, maximum component-scaled error
**1.639569408872938e-15**, fixed acceptance bound `3e-12`.
Each case has 4x4x4 periodic cells and 240 faces across three axes, comparing
mass, three momenta, total energy and seven species fluxes. The three fixtures
are uniform, multidirectional pressure/velocity/temperature variation, and
additional independent mixture variation. Each uses dt=0 and dt=1e-5 with
spacing 0.25. Scaling is the maximum absolute flux for each component across
all faces in the case (reference and actual), floored at one.

Fortran stream mode calls production `compute_reactive_ppm_ctu_fluxes_3d`
with `pelec`, zero viscosity, no flattening and no contact steepening. The
reference runs the full `pc_umeth_3D` periodic zero-source pipeline plus final
species normalization. Pressure-divergence and internal-energy flux are computed
by the reference but are not observables in the comparison.

The original normal/CTU unit passed **1/1**, 6.26 s test / 6.31 s total.
Provenance, Python syntax, project contract and whitespace passed. No failed or
skipped checks. GNU Fortran/C++ 13.3.0, Python 3.12.3, existing Ninja Debug
tests-enabled `build/spray-amr-reactive-debug`. Build warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_normal_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tools/check_pinned_ppm_ctu.py --reference /mnt/sn850/Projects/gtt2/external/PeleC --executable build/spray-amr-reactive-debug/tests/test_reactive_ppm_normal_3d --work-dir build/ppm-pinned-ctu-checks/oracle
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_normal_3d$' --output-on-failure
python3 -m py_compile tools/check_pinned_ppm_ctu.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/ppm-pinned-ctu-checks/{build.log,comparison.log,ctest.log,unit.log}`
and `oracle/` with pinned sources, generated driver, compiler log, inputs,
thermo data, outputs and report. The report records compiler command/version,
source/shim/EOS-adapter/executable hashes. Reference:
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

This establishes bounded periodic zero-source assembly agreement after the
[independent-species projection](ppm-mixture-projection-20260926.md) and
[passive upwind/normalization](ppm-passive-flux-20260926.md) repairs.
Artificial viscosity, shock flattening, physical boundaries, forcing, EB,
selected mechanisms, MPI/application execution, installed audit and physical
validation are outside this oracle. Complete PeleC parity remains unproven.
