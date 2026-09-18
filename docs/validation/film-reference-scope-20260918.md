# Stationary film reference scope — 2026-09-18

Repository base: `3ad36c5`. Pinned PelePhysics reference:
`6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`.

## Finding

Lateral wall-film transport is a future model extension, not an unported
operator in the pinned spray update path. This classification does not remove
the extension from future work or establish complete spray/film parity.

In [SprayParticles.cpp](https://github.com/Pele-Suite/PelePhysics/blob/6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b/Source/Spray/SprayParticles.cpp#L362),
positive film height sets `is_film`; the film branch calls `calculateFilmSource`.
The position update at lines 449–452 explicitly requires `!is_film`.
In [WallFilm.H](https://github.com/Pele-Suite/PelePhysics/blob/6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b/Source/Spray/BreakupSplash/WallFilm.H),
the source routine changes temperature, component inventory and height while
keeping diameter fixed. It supplies no tangential-position update. Its comment
also identifies the fixed-diameter normal-loss assumption as a limitation.

PeleF's `spray_wall_film` retains position/normal and stationary thermodynamic
state; its existing wall-film evolution follows this stationary scope. Adding
lateral momentum/transport would need new equations, state and qualification.
No invented lateral closure was added as part of this reference port.
Existing differences in wet-impact closure, mixture thermodynamics, transfer
integration and physical qualification remain separate open questions.

## Provenance and checks

Both files were freshly retrieved from the pinned raw GitHub URLs. SHA256:

- `Source/Spray/SprayParticles.cpp`:
  `0e376f77df434ae7eec9fba39ab67af5ad4c20fb4c46529a9ff35b6d978b4d05`
- `Source/Spray/BreakupSplash/WallFilm.H`:
  `6b8f1156555058059ab21c023444a662c3bc1f7d23f3173b83ccc92eeb0e51e5`

The film hash matches the existing manifest. `references/spray_wall_impact.json`
now also pins the position-update source. Retrieved files and URL/hash records
are in `build/film-reference-scope/`. Mechanism provenance, manifest JSON parse,
Python syntax, project contract and whitespace checks passed. No solver build,
numerical tests, physical validation or full completion audit was performed for
this documentation/provenance change.
