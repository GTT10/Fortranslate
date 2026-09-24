# Shared gas-species mapping for liquid components — 2026-09-25

Implementation `aebdb94a79db94d4c975070f129aba9f5a9692e0`, base `dfd5980`.
Pinned PelePhysics `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`, `Drag.H`:
SHA256 `bcc2bcb2f01b8cedfeb92d0bfc9ae2c8e723e77285f7011aeb77d5b138a29345`.
The source bytes are retained in `build/spray-property-reference-audit/Drag.H`
and listed in `references/spray_property_scope.json`.

## Numerical change

Components may now share a gas species if their gas names agree. Liquid records
and component inventories remain separate. Surface contributions are summed
before skin composition is formed, active flags are combined, and parcel and
wall-film evaporation sources accumulate into shared gas species instead of
overwriting prior contributions. Existing inventory/boundary accumulation was
already additive.

The reference at `Drag.H:103–116` first computes ambient allocation
`Xgn=Xdn*Xgi`, dividing by the sum of liquid mole fractions in that gas-species
group when the group contains another nonzero contribution. It compares
**pure** psat with `Xgn*P`, then retains `Xdn*psat/P`. PeleF now applies this gate
to its Raoult surface fractions. For a unique mapping this also differs from
the previous `Xdn*psat/P > Xgi` test. Shared vapor is then collapsed as in the
reference's subsequent gas-species accumulation.

Two old expected behaviors changed intentionally. In the wet-mixture fixture,
both components pass the reference gate, but their combined surface mass
fraction is below the ambient fraction, so total evaporation is zero. In the
trace-film fixture, the allocated ambient pressure is tiny and the trace stays
active; existing trace depletion releases its full residual mass. Updated tests
assert both the gate outcome and the resulting source, without changing any
conservation tolerance.

Mixture evolution identity is now `mixed-liquid-components-v3-vapor-mapping`.
Regular serial/MPI, MPI AMR and serial/MPI EB application contexts bind it.
Old model contexts are intentionally incompatible. The gas-species molecular
weight still supplies each mapped liquid component's molecular weight; arbitrary
independent liquid molecular weights, full reference skin/transfer parity and
experimental validation are not established.

## Completed focused checks

GNU Fortran 13.3.0, Debug/Ninja, Python 3.12.3; unchanged fixed tests-enabled and
selected H2/O2 tests-disabled builds. Fixed core and regular application built.
Selected regular, MPI AMR and MPI EB applications built. Compiler warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_core pelef_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_spray_3d_selected pelef_mpi_amr_spray_3d_selected pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_(many_to_one|mixture_rates|mixture_wall_rates|mixture_exchange|mixture_wall_exchange)$' --output-on-failure --output-log build/shared-vapor-checks/ctest-final.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/spray/test_spray_application.py --exe build/spray-amr-reactive-debug/pelef_spray_3d SprayApplication.test_many_to_one_restart
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/spray/test_spray_application.py --exe build/spray-amr-reactive-selected-debug/pelef_spray_3d_selected SprayApplication.test_many_to_one_restart
python3 -m py_compile tests/spray/test_spray_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final core checks **5/5 passed**, total **0.06 s**. New coverage includes shared
dry/humid vapor equivalence, inactive-component ordering, additive parcel/film
gas sources and retained liquid component masses. Initial core run failed the
two obsolete humid/trace expectations described above; logs are retained and
the test target was rebuilt before the final run. No tolerances were loosened.

Fixed/selected application checks each passed **1/1**, in **0.473/0.481 s**.
They use two components mapped to H2O with different liquid properties, check
existing mass/energy/momentum bounds, exact gas/parcel CSV restart equality and
changed-property rejection. No skips. Provenance, Python syntax, project contract
and whitespace passed. Artifacts: `build/shared-vapor-checks/`, including build,
initial/final CTest, copied LastTest and fixed/selected application logs.

No distributed execution of shared-vapor liquids, full matrix, installed-binary
audit, external solver run or physical mixture validation was performed here.
