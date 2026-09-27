# Serial inviscid rotating MOL time integration — 2026-09-27

Source `d83162d42121b2aa0d3aa62fd216bcb48d616dfa`, main, base `ccb749ff`.

## Change

The optional frame now reaches common-step and recursive serial MOL SDC.
Both old and predicted endpoint momentum states feed centrifugal/Coriolis
sources, with trapezoidal integration into the existing source ledger.
Rotation adds no mass, species or rotating-total-energy source. A separate
applied force retains its existing work term.

Post-update EOS, EB pressure ledgers, peer/reflux corrections and restriction
use receiving-cell rotating potentials. Existing private candidates retain
transactional state/ledger publication. Nonrotating call sites remain valid.

This stage accepts **nonreacting, inviscid MOL** with a frame. CTU, reaction or
any transport flag with a frame is explicitly rejected until those operators
carry the energy convention. This is an implementation boundary, not a reduced
completion goal. Applications do not yet expose frame input or restart identity.

## Completed verification

GNU Fortran 13.3.0 Debug, fixed core and isolated heptaneSimple selected serial
application builds. Existing build has MPI/spray ON and CTest targets/SUNDIALS
OFF. Manual tests use `-g -fcheck=all -ffree-line-length-none`. Build/link commands
unset both library-path environment variables. Existing executable-stack
linker warnings remain.

Extended `test_rotating_eb_tree_3d` retains the earlier
[transfer checks and fixture](rotating-eb-tree-transfer-20260927.md): two levels,
three all-fluid 4³ patches; pure N2, 300 kPa, density 1 kg/m³, velocity
`(1,2,-1) m/s`; z rotation at 100 rad/s through `(10,10,0) m`, with negative
rotating total energy. Each new case advances one root step of `1e-9 s`:

| Mode | Flux | Redistribution | Maximum composite relative residual |
| --- | --- | --- | --- |
| Common | acoustic | none | `3.5568237000791222e-15` |
| Common | LF | flux | `3.8634042907883517e-15` |
| Recursive | acoustic | none | `3.5568331778576566e-15` |
| Recursive | LF | state | `3.8634107237124962e-15` |

Residual is `max(abs(U_after-U_before-boundary-wall-chemistry-source) /
max(1,abs(U_before)))`, from composite integrals, gate `3e-14`. All four passed;
recursive fine levels take two substeps. Source mass/species/energy were exactly
zero and the centrifugal momentum sign was correct. The two common cases also
compare the source ledger with independently evaluated z-axis formulas using
the old state and Euler predictor; absolute error passed `1e-15`.

The combined executable passed 6,207 assertions, mostly the existing per-cell
transfer checks; this count is not a physical-case count. Existing
`test_serial_sdc_species_3d` passed its four ordinary CTU/MOL acceptance/rollback
cases. Provenance, artifact Python syntax, project contract and diff whitespace
checks passed. No repository Python changed. No test failed in this increment.

## Artifacts and remaining work

Exact commands and logs: `build/rotating-sdc/commands.txt` and neighboring files.
Manual executable SHA-256:

- `check`: `98048af3417bc571d4714976a44da521050658ea8cc31da1cf3114190a1d3eb1`
- `default`: `1b3c1d54ffa19613e1e50bb67ea42fd3724f502a06a840905e0f01e3e57b080b`

Selected serial application SHA-256: `ad299245995bfac5b406ab6f14dfc488e9c5928c747dafb31fbd15af32aa9d0b`.
It was built, not run. No full CTest matrix,
installed-binary audit, new compiled PeleC comparison, rotating cut-cell
trajectory or physical validation was performed. Reaction/transport integration,
MPI owner-local stepping, conservative remapping, timestep selection and
application/checkpoint frame identity remain open. The earlier 1000 rad/s LF
absolute patch-conservation failure remains open. The frozen long Jet run was
observed live and preserved.
