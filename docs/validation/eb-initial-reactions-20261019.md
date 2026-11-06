# Fresh EB MOL reaction-history initialization, 2026-10-19

Gas EB tree applications now initialize their first reacting MOL predictor
history independently. This extends the [regular tree initialization](mol-initial-reactions-20261019.md)
to embedded boundaries without advancing the initial gas state.

## Numerical change

The serial and owner-local MPI initializers reuse the existing EB reaction
stage with private gas, temperature, internal-energy and history copies and
zero spatial forcing. Only fluid cells enter the reactor. The complete new
history is published only after every level and rank accepts it; failed
initialization preserves the caller's history. Solid-cell history is zero.

The root uses the first attempted coarse interval; recursive finer levels
use successively divided intervals, while common clocks keep the root
interval. Ordering follows the pinned `PeleC::post_init`: root reaction,
volume-weighted restriction of still-zero fine history, then fine reaction.
No new average-down of public gas or temperature is introduced.

The application connects initialization only to fresh reacting gas MOL.
It runs after timestep/event limits and uses the selected native/CVODE
backend. A rejected first gas step restores the previous history before
retrying at a shortened interval. The first accepted step reports one
initialization; checkpoint continuation uses its saved history without
reseeding. Existing checkpoint layouts and contexts are unchanged. The
separate particle startup is outside this change.

## Completed integration checks

- One three-rank unit CTest passed in 0.48 s. Common and recursive clocks
  initialize all 240 fluid cells once, including cut cells; invalid gas in
  solid cells is never passed to the reactor. Serial/MPI histories match
  exactly, an idle rank performs no reaction work, public gas/T/UEINT stay
  unchanged, and a late fine-level failure restores the entire history.
  The synthetic H2-to-H reaction checks the interface and rollback, not
  physical chemistry.
- The existing hot CVODE EB-AMR application checker passed with recursive
  clocks. Five successful calculations cover serial continuous, reaction-off,
  two-rank continuous with rebalancing, two-rank checkpoint stop and three-rank
  continuation. Two additional invalid output policies were correctly refused
  before integration. No unexpected failure or skip occurred in these checks.
- Final and intermediate fields are byte-identical between serial, two-rank
  continuous and changed-rank restart at their common physical times. Each
  fresh reacting run reports one initialization; reaction-off and restart
  report none. Independent CSV inventories, species closure, H/O/N chemistry
  balances and NASA7 UEINT audits pass. Maximum scaled budget residual is
  `6.941e-16` (the existing check uses a scale floor of one); hydrogen consumed
  is `2.29522694e-4` of the initial inventory, or 0.02295%, over 300 ns.

The application evidence is
`app-runs/reacting_restart_1vv7zsj7/checks.json`; exact commands and timings
are in `commands.json`. Its seven processes took 36.72 s including checker
work. This is the selected H2O2 mechanism case; LiDryer is used in the unit
interface fixture and the separate fixed-reference comparison.

## Fixed-reference check

The unchanged pinned PeleC EB-C3 and the new selected-LiDryer Fortran binary
were both rerun for the existing uniform 2000 K planar-EB case. Four times
(0, 100, 200 and 300 ns) pass all existing field/geometry/inventory gates
without tolerance changes. Maximum temperature difference is `7.7761e-8 K`,
species mass-fraction difference `4.8163e-10`, and hydrogen-consumption relative
difference `2.2235e-5` (0.002224%). The reference consumed `4.478803865e-5` of
its initial hydrogen; Fortran consumed `4.478704281e-5`.

The result is `reference-runs/matched_eb_cfe_pmak/checks.json`, including
reference input/source contract, binary hashes, commands and physical times.
This independently rerun uniform single-level case does not establish
spatial-gradient EB-AMR parity. Its checker took 22.14 s: reference integration
2.72 s, Fortran 18.68 s and four plotfile exports. The entire focused unit,
application and reference verification took about one minute after compilation.

Both LiDryer and H2O2 mechanism provenance/regeneration preflights passed.
Python syntax for the application and reference checkers, project contract
and `git diff --check` passed. No full matrix or unrelated regression rerun was
performed. The four selected serial/MPI gas frontend targets and the new unit
target built successfully; fixed-mechanism executables were not rebuilt in
this increment.

## Identity and scope

Branch `main`, base `0a0e00ebe9138145fbb0e0e745665d0d768df009`, with the existing
dirty working tree preserved. GNU Fortran 13.3.0 Debug, Open MPI 4.1.6,
CMake 3.28.3 and SUNDIALS 7.2.0; compiler/MPI/backend configuration unchanged.
Artifacts are under `build/eb-initial-reactions-20261019/`: `before/`,
`increment.patch`, `source/`, `source-identity.json`, `CMakeCache.txt`,
`commands.json`, logs and unique application/reference result directories.
The source snapshot records 958 source/input files including the pre-existing
uncommitted numerical implementation. `increment.patch` isolates this change
from that earlier work; added files are retained in `source/`.

Fresh reacting MOL results can change from the former zero-history startup.
This increment is not a full configuration matrix, installed-binary audit,
spatial-gradient EB parity, mesh/time convergence or physical ignition
validation. Rotating, particle and LES startup qualification remain outside
the executed cases. The nonreacting MOL2 conservation-policy difference is
unchanged.
