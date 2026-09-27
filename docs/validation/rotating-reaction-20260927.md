# Rotating-energy SDC chemistry — 2026-09-27

Source `59f0d1ca344a9ed9766ea79c42aec1a38076446e`, main, base `0e366625`.

## Change

The 3D forced SDC reaction adapter accepts a potential field. It recovers old
and nonreacting endpoint internal energies by adding the corresponding `rho*Phi`
to rotating total energy before reaction integration. After reaction, it uses
the resulting density to restore rotating total energy. The inferred energy
reaction source follows the same convention. This includes nonreacting density
forcing; treating `rho*Phi` as a constant would give the wrong thermal forcing.

Serial EB tree chemistry supplies each receiving cell's potential. Common and
recursive rotating MOL SDC now accept chemistry. Diffusive transport remains
explicitly rejected with a frame until its energy convention is connected.
No reaction mechanism, reactor tolerance or generated source changed.

## Completed checks

GNU Fortran 13.3.0 Debug fixed core; isolated heptaneSimple selected serial build.
The existing configuration enables MPI/spray, disables CTest targets/SUNDIALS.
The selected executable was built, not run. Manual links use `-g -fcheck=all
-ffree-line-length-none` with both library-path environment variables unset.

- New `test_rotating_reaction_3d` passed. A 1 kg/m³, 600 kPa mixture with mass
  fractions H2=.028, O2=.222, OH=.0001, N2=.7499 and velocity `(10,-20,30) m/s`
  advances `1e-7 s` with potential `2e6 J/kg`, negative rotating total energy,
  N2/density source `3 kg/m³/s`, momentum source `(2000,-1000,300)` and rotating
  energy source `3e5 J/m³/s`. Ordinary comparison adds `Phi*S_rho` to its
  energy forcing. After converting back, states agreed exactly (gate `3e-11`);
  temperature passed the same relative gate. H2 density changed by
  `-2.0927391772795206e-5 kg/m³`, demonstrating actual reaction in this cell test.
- Zero potential reproduced ordinary state, temperature and inferred reaction
  source exactly. Invalid potential retained all three output arrays exactly.
- Extended tree test passed six modes: the previous four inviscid nonreacting
  modes plus two-iteration reacting common/acoustic and recursive/LF modes.
  Those tree cases retain **pure N2** and verify the coupled chemistry path,
  not chemical conversion. Their composite relative residuals were
  `3.6902944693328716e-15` and `3.9084276956431077e-15`, below unchanged `3e-14`.
  The executable reports 6,217 assertions including prior per-cell transfer checks.
- Existing ordinary `test_serial_sdc_species_3d` passed four CTU/MOL
  acceptance/rollback cases. Fixed/selected builds, provenance, artifact Python
  syntax, project contract and diff whitespace checks passed.

Initial manual compilation from the repository root twice encountered an old
root-level `reactive_3d_mod.mod` (dated September 7) instead of the current build
module. Compiling from the artifact directory resolved this; old files were
preserved. An initial unseeded H2/O2/N2 fixture failed the explicit chemical-change
assertion at this short interval. The committed fixture includes OH so the test
actually exercises chemistry; solver/test tolerances were not loosened. Existing
executable-stack linker warnings remain.

## Artifacts and remaining work

Commands and logs: `build/rotating-reaction/commands.txt` and neighboring files.

| Executable | SHA-256 |
| --- | --- |
| `check` | `612e438c10916a3a85ab081351d42e26abd71ce331b9c234e026c4468b9e4624` |
| `tree` | `87d44d8a21346b25f04fbd3f927774e5e51c9fdfa717fdb2fe4c14481df43bcd` |
| `default` | `374d92222d39cc4e7af1b5844b625d5156c65ebed6773652535b129dc2120fbe` |

No full matrix, installed-binary audit, compiled external comparison, physical
validation or reactive rotating-flow trajectory was run. Transport, MPI,
remapping, timestep selection and application/restart frame identity remain
open, as does the prior 1000 rad/s LF absolute-conservation test condition.
The frozen long Jet continuation was observed live and preserved.

Selected serial executable SHA-256: `c017097f7db8764a0894608821bceeb560dfc02fd4e1214fa79c631b903f8c11`.
