# EB spray periodic field exchange and coupling — 2026-09-25

## Source and behavior

Implementation `1264d55`, based on `705d54d`, local main.
The EB cloud-in-cell stencil, sampler and conservative deposition accept an
optional three-axis periodic mask. Stencil distances and face connectivity use
unwrapped donor indices; storage indices wrap into the physical mesh. This
preserves local weights and connectivity at seams, including corners and
repeated donors on small domains. Covered cells remain excluded. Opposing
geometry faces and nodal level-set endpoints must agree before exchange.

The EB particle interval forwards its existing outer-face periodic flags to
both sampling and deposition. The allocated serial reacting-spray wrapper now
accepts paired periodic gas faces and forwards those same faces to particle
drift. Inflow/outflow mapping is unchanged. Explicit periodic faces with the
impact/film model are rejected; that model currently supports only outer
outflow. The application still restricts explicit input to inflow/outflow.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug, Python 3.12.3; existing fixed tests-enabled
build/spray-amr-reactive-debug. No toolchain/backend changes.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_eb_field test_spray_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_(field|reactive|reactive_ctu)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**3/3 passed, 0 failed, 0 skipped**, 2.96 s total: field 0.01 s,
legacy coupling 0.30 s, CTU coupling 2.61 s. No failed test attempts.

The field fixture samples a varying y/z field at both equivalent periodic
corners beside the x-normal cut plane. Four donors have equal weights and
produce the analytic value 27.5 from both sides. Deposition reaches both seam
sides, preserves integrated inventory to the existing 1e-14 tolerance, and
leaves covered storage unchanged. Incompatible periodic geometry and nodal
seams reject without publishing samples.

The coupled fixture passes a moving parcel through the periodic y face with
heat/evaporation, chemistry and transport active. The parcel retains positive
mass and appears on the opposite side. Gas-plus-particle inventory closes with
the existing 3e-12 relative tolerance. Existing open-boundary, force and failure
rollback cases also pass. No tolerances or golden data changed.

Provenance, Python syntax, project contract and diff checks passed. Mechanism
SHA remains 0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4.
Artifacts: build/eb-spray-periodic-field-checks/build.log and LastTest.log.
Existing real-comparison compiler warnings remain.

## Limits

These are fixed-mechanism serial library checks. Periodic application input
and restart, selected runtime, MPI/tree dispatch and periodic impact/film
coupling remain outside this evidence. No full configuration matrix, installed
binary audit, external solver parity or physical spray validation was run.
