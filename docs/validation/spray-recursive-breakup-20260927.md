# Recursive spray breakup and computational splitting — 2026-09-27

Implementation `67dbdc39907c03844bf4572cc18a0aaf88c39b3a`, based on `ef76b767`,
local main, version 0.246.0.

## Change and reference boundary

The recursive particle halo kick now calls the existing ETAB/KHRT diameter
update after frozen-gas heat/mass/momentum exchange and before drift. It uses
mixture density, heat capacity and energy plus the mixture-film gas properties,
as the existing common-step coupling does. Component masses scale with the
updated droplet mass; the breakup kernel preserves represented liquid mass
through the corresponding multiplicity change. Physical, virtual and ghost
populations all update their diameter/history. Intermediate SDC iterations
still perform no particle update.

After physical collection at outer synchronization, optional
`maximum_multiplicity` splitting uses the existing MPI allocator for globally
unique IDs and component inventories. Its outputs remain private until the
entire gas/particle transaction succeeds. `split_created` reports the new IDs;
application diagnostics now receive that result. Ranks agree on the split
threshold before branching. KHRT shedding children remain explicitly rejected
by the recursive wrapper, halo kick and application; films remain unsupported.

Pinned PelePhysics `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`,
`Source/Spray/SprayParticles.cpp`, updates breakup following particle exchange
and before motion, including copies. Its creation of physical children is gated
by active/moving particles. This change reuses PeleF's existing diameter kernels
and final computational splitting. It does **not** reproduce upstream's complete
child-creation timing, tangential child velocities or KHRT shedding lifecycle.
Those remain separate integration work; this is not exact upstream spray parity.

## Independently completed checks

Artifacts: `build/spray-recursive-breakup-checks/`.
Existing fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`,
GNU Fortran 13.3.0, `/usr/bin/mpifort` and `/usr/bin/mpiexec` OpenMPI 4.1.6.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d test_spray_mpi_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc_recursive_(breakup|khrt)_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final tests **2/2 passed**, no failures or skips: ETAB 13.36 s, KHRT 13.32 s,
26.73 s total. Maximum composite scaled residuals were respectively
`7.1299975393456751e-16` and `1.1112755134928231e-15`, against the unchanged
`3e-12` criterion. Provenance, final build, Python syntax, project contract and
whitespace checks passed. An initial test build exceeded the Fortran line-length
limit; the two assertions were wrapped before testing. Its log is preserved in
`build-line-length-failed.log`. Existing executable-stack linker warnings remain.

Each two-rank case advances a two-level reacting/transporting EB gas with a
mixed-liquid physical parcel and recursive copies. It checks early gas rejection,
late split-capacity failure, unchanged gas/temperature/particle inventories and
IDs on failure, then a successful retry. Accepted output respects the split
multiplicity limit, has the expected final owner and new-ID count, and closes
the coupled mass/species/momentum/energy inventory. The ETAB fixture starts with
deformation 1.1 and explicitly requires droplet mass reduction beyond evaporation;
the KHRT fixture exercises its update path and computational splitting without
claiming an independently verified KH/RT instability event.

No accelerated-breakup combination, application runtime, selected-mechanism
matrix, rank-changing restart, three-level breakup or physical validation was
run in this focused check. Application wiring was compiled.
