# Rotating-energy serial tree transfer — 2026-09-27

## Implementation

Source `7076447b3fe91ec007c15a60d9ba13a19c744853`, main, base `5df0f6bc`.
The optional rotating frame now reaches serial EB tree initialization,
restriction, halo filling and MOL spatial rates. Underlying PCM prolongation,
restriction, patch halos and peer overlays accept the same frame.

Conserved transfer averages or copies rotating total energy without substituting
ordinary total energy. Thermodynamic recovery uses the **receiving cell center**.
Thus fine-cell temperatures after PCM transfer can differ from the coarse
cell temperature while the conserved inventory remains unchanged. The existing
transactional publication of states and temperatures is preserved.

Prescribed and Mach inlet APIs retain ordinary relative total energy at their
interface; the halo adapter converts to rotating energy at each receiving cell.
Mach reconstruction temporarily converts its interior argument to ordinary
energy. No rotating source or time integration is enabled by these APIs alone.

## Completed checks

GNU Fortran 13.3.0, Debug fixed core and isolated selected serial build, with MPI
and spray enabled, CTest targets and SUNDIALS disabled in the existing build.
The selected build uses heptaneSimple. Manual tests link with `-g -fcheck=all
-ffree-line-length-none`; build/link environments unset `LD_LIBRARY_PATH` and
`LIBRARY_PATH`. Exact commands are in `build/rotating-tree/commands.txt`.

- New `test_rotating_eb_tree_3d` passed: a 4³ root, two adjacent 4³ children,
  two levels, all-fluid geometry in a 1 m cube. Pure N2 starts at 1 kg/m³,
  300 kPa and relative velocity `(1,2,-1) m/s`. Rotation is 100 rad/s about z
  through `(10,10,0) m`, making every root total-energy density negative.
- PCM initialization followed by restriction reproduced root conserved values
  exactly; recovered root temperature agreed within `1e-12` relative.
- Each of three patches filled a 10³ halo volume. Independent EOS calls using
  each receiving position reproduced the returned temperatures exactly.
  Explicit peer and parent-donor conserved-state comparisons also passed.
- Acoustic/LF × order 1/2 MOL tree spatial calls succeeded with integration
  weight `dt=1e-9 s`. These are four operator calls, not evolved trajectories.
- Prescribed and Mach inlets passed with one and three halo layers. Prescribed
  energy subtraction was exact; the Mach inlet recovered 900 K, 300 kPa and
  Mach 0.2 within `1e-12`. Invalid-frame rate/flux/restriction rollback passed.
- The executable reports 6,177 assertions, mostly per-cell input/EOS checks;
  this is not 6,177 independent physical cases.
- Existing `test_eb_reactive_transfer_3d` and `test_serial_sdc_species_3d`
  passed (the latter covers four CTU/MOL negative/positive species-source cases).
- Fixed core and selected serial application builds, provenance, artifact
  Python syntax, project contract and `git diff --check` passed. The selected
  application was not run. No repository Python changed.

An initial build failed on two extended Fortran lines exceeding 132 columns;
continuations corrected the source and the completed builds passed. Existing
executable-stack linker warnings remain. No tolerance or golden change was made.

## Artifacts and limits

Logs and manual executables: `build/rotating-tree/`.

| Executable | SHA-256 |
| --- | --- |
| `check` | `7b8807526a5012bc0c4bcdc2c5f8bb564995375eaa16b9edf05f36dc7b9845ad` |
| `transfer` | `d782feaa52823fee963115ff63aaa3095e502949bd313fba6a4244ad4672e184` |
| `default` | `43ba8b3bbcf91391b8deb3a8da2b33d33fc0c38c0747a9b983bdd92a5f5da3da` |
| selected serial | `d01e1d0b4d9453f44fc5f1091744bb7ea4e07ef630d85169b7d67385091a5c1e` |

This is serial transfer and spatial-operator evidence. MPI owner-local paths,
old/new rotational sources, reaction/transport, time integration, reflux/remap
recovery and application/checkpoint frame identity still need connection.
No new compiled external-solver comparison, rotating cut-cell tree trajectory,
full CTest matrix or physical validation was performed. The prior
[1000 rad/s patch conservation-gate failure](rotating-mol-patch-20260927.md)
remains unresolved. The frozen long Jet continuation was observed running and
was not modified.
