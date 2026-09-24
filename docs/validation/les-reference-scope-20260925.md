# LES reference scope audit — 2026-09-25

Audit base: `daead4c`. Fresh source downloads use the frozen PeleC commit
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`; they do not follow a moving branch.
This audit corrects the classification of remaining work, not the solver.

## Embedded boundaries

[Source/LES.cpp](https://github.com/AMReX-Combustion/PeleC/blob/bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8/Source/LES.cpp)
checks the cell-box `FabType` in all four closures. The Smagorinsky (line 270),
Dynamic Smagorinsky (414), WALE (658) and Vreman (765) paths invoke
`amrex::Error` when that type is not `regular`. A subsequent covered-box check
does not establish cut-cell support; the preceding non-regular guard applies.
The Dynamic routine begins at line 341.

Therefore LES closure on non-regular EB cells is **not an omitted pinned
PeleC operator**. It remains a possible future model extension requiring its
own cut-cell filtering, stress/heat flux and qualification design. Earlier
PeleF status wording that listed Dynamic EB alongside missing reference
functionality was too broad.

An EB factory containing regular cells is not the same as a cut-cell LES
operator. PeleF's regular-grid serial/MPI AMR LES path is the implemented
alternative for LES. EB hydro, molecular transport, chemistry and spray
capabilities remain separate; this finding neither removes their requirements
nor claims that existing EB applications expose LES selection.

## Time integration remains required

[Source/Advance.cpp](https://github.com/AMReX-Combustion/PeleC/blob/bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8/Source/Advance.cpp)
accepts the AMR iteration/cycle indices at each level and dispatches source
evaluation during the advance. [Source/Sources.cpp](https://github.com/AMReX-Combustion/PeleC/blob/bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8/Source/Sources.cpp)
routes `les_src` to old/new LES construction (lines 25–26 and 82–83).
In LES.cpp, old/new source construction uses time and dt; the new-source
reflux weight depends on the sub-iteration, and `updateFluxRegistersLES`
accumulates `reflux_factor * dt` (lines 217–229).

This establishes a level/source/flux-register integration responsibility, not
numerical parity between PeleC's integrators and PeleF SSPRK2 splitting.
PeleF's current LES operators use a common interval across the tree. Separate
recursive hydro/transport support does not prove recursive LES completion.
The next implementation work remains per-level LES stepping with temporal
coarse donors and interval-weighted reflux, followed by MPI/coupled dispatch.

## Reproducibility

Downloaded files reside in `build/les-eb-reference-audit/`:

| Source | SHA-256 |
| --- | --- |
| LES.cpp | `9106d8bb7623d6ad1200b0bc4dab8020d62b6d4ce557c557f18db09f175c87f1` |
| Advance.cpp | `acde092d71d158da8cafc87cc44a2e06b124af82f472e4c7050cc47e395ecbc0` |
| Sources.cpp | `14489cdef3e6b52d73346de7d88839f02431a72628d9d0235f7b1060775f5bf0` |

The fresh LES.cpp bytes match `build/les-reference/LES.cpp` exactly. Files were
retrieved using Python `urllib.request.urlopen` from the corresponding
`https://raw.githubusercontent.com/AMReX-Combustion/PeleC/<commit>/Source/<file>`
URLs and hashed with `hashlib.sha256`.
Mechanism provenance, project contract, Python syntax of the contract checker
and whitespace checks passed. No numerical calculation, build, install audit
or physical validation was performed for this documentation-only audit.
