# Spray property scope audit — 2026-09-25

PeleF source reviewed: `03d6097`. Frozen PelePhysics:
`6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`, as recorded by the PeleC baseline.
This is a source/documentation audit, not a numerical implementation or parity run.

## Reference boundary

`SprayPropertiesGCM.H:310–315` rejects `PELE_SPRAY_GCM` unless
`PELELM_USE_SPRAY` is defined and explicitly identifies PeleLMeX as the supported
application. GCM is therefore not a missing property path for this pinned
PeleC port. The existing PeleMP-style property path is the alternative.

`SprayPropertiesPeleMP.H:94–103` defines the Watson correlation
`L(T) = Lref * ((Tc - T)/(Tc - Tref))**0.38`. Its comment mentions manifold
EOS, but that does not delimit all uses: `SprayFuelData.H:135–155` evaluates
it at the normal boiling temperature and uses that latent heat in the
pressure-dependent Clausius boiling estimate. `SprayParticles.cpp:387` calls
this routine before the film/ordinary-droplet branch, without a manifold guard.

The ordinary `calcVaporState` in `Drag.H:50–114` instead evaluates gas enthalpy
at the input droplet temperature and calls `latent_i` with the capped surface
temperature. That result feeds the saturation pressure. The manifold
specialization at line 172 uses the correlation. A blanket Watson substitution
would therefore change the ordinary EOS path incorrectly.

## Concrete remaining port work

PeleF's `spray_parcel_mod.F90::saturation_mole_fraction` currently uses one
`reference_temperature` both for stored liquid energy and for the saturation
reference state. Its optional `surface-cap` uses the supplied constant latent
heat. The frozen PeleMP source distinguishes a latent-property reference
temperature from the normal boiling temperature. Matching that path requires:

- Separate, validated boiling and thermodynamic reference inputs, preserving
  the current defaults and binding any new settings into restart context.
- Watson evaluation at the normal boiling point for the pressure-dependent cap.
- A separately specified enthalpy-derived vapor-state latent heat for ordinary
  EOS; do not silently change the conservative stored-liquid-energy convention.
- Focused independent formula checks and a coupled continuation/conservation
  check once the new path reaches application input.

This audit does not establish physical mixture boiling, supercritical behavior,
or reference time-integration parity. The older surface-cap evidence remains
valid for its explicitly constant-latent scope.

## Documentation reconciliation

The roadmap and current-status table no longer list all MPI/AMR particles,
distributed feedback, film ownership or application restart as unimplemented.
The existing AMR film and MPI EB regrid evidence linked there documents bounded
implementations. Those historical test counts were not rerun here. Broader
qualification remains open; lateral film transport remains an extension beyond
the pinned stationary-film model.

## Provenance and verification

All five reference files were extracted with
`git -C /mnt/sn850/Projects/gtt2/external/PeleC/Submodules/PelePhysics show
6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b:Source/Spray/FILE`.
The checkout's current files were not used as authority. SHA256 values are in
[`spray_property_scope.json`](../../references/spray_property_scope.json);
extracted files and check logs are in `build/spray-property-reference-audit/`.

Verification commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Python 3.12.3 was used. Manifest parsing/hash comparison (5/5 files) and the
listed checks passed, with no failures or skips. Check output is retained in
`build/spray-property-reference-audit/checks.json`. No solver builds,
numerical tests, full configuration matrix, external solver comparison or
physical validation were run for this documentation-only change.
