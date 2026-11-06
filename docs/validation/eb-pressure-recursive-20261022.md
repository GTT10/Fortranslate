# Recursive pressure clocks and native MOL controls, 2026-10-22

The recursive candidate prepares covered-parent pressure only at the start
of each synchronized subtree clock. It evaluates fluid EOS pressures and
restricts them with fluid-volume weights within that subtree, then reconstructs
only the current level's covered gas/held/temperature representation on the
private candidate. Exposed valid leaf cells are untouched by this operation.
The initialization is before clock-history capture; pressure differences at
floating-point roundoff are not reapplied when entering a saved clock.

During provisional level updates, spatial halos use that level's current EOS
pressure without restriction from descendants at other times. For ancestor
donors, pressure is interpolated from the two parent endpoint EOS pressures,
using the same ancestor fraction and stencil as the material reconstruction.
The material interpolation remains conservative, followed by the existing
auxiliary pressure reconstruction. MPI owns and exchanges the pressure field
with the same donor schedule and collective admission as gas/held data.

This leaves the common-clock policy unchanged. The pressure policy remains
opt-in, requires both momentum options, held MOL with no redistribution, and
a nonrotating frame. Ordinary gas application subcycling now reaches this
candidate; saved policy identity and checkpoint envelopes still bind it.

## Verification boundary

Artifacts: `build/eb-pressure-recursive-20261022/`.
The temporal kernel test distinguishes endpoint-pressure interpolation from
EOS of a mixed nonisothermal state. Application checks exercise fixed/selected
recursive restart, transport, chemistry, regrid and changed ranks. Dedicated
pressure replay/retry tests exercise the existing forced-refusal fixtures.
Completed results, source hashes and commands are recorded with the artifacts.

A new `--audit-native` probe mode uses the production common/recursive held-MOL
integrator on the saved smooth-isothermal and imported nonisothermal controls,
with chemistry, molecular transport and forcing disabled. It compares
velocity-only and paired-pressure policies over identical imported states.
Each step writes independent composite inventories and returned physical and
numerical-energy changes, and checks all component budgets. Pressure error
is observed, not reset on valid leaf cells. This replaces private SSPRK2
integration as evidence for the behavior of the ordinary numerical core.

No physical convergence or fixed-PeleC qualification follows from an
application/replay test. In particular, eliminating nonisothermal pressure
error by forcing valid-leaf pressure would change the energy method and
requires an explicit numerical-energy ledger. No such forcing is introduced.

## Native production-core observations

The frozen probe in `native/frozen-probe` runs four 50 ns root steps, two MOL
iterations, the saved 8-cell-per-axis EB fixture, and no chemistry/transport.
Both policies enable separate momentum and consistent velocity. Pressure
is measured on exposed fluid leaves relative to the initial 101325 Pa.
Values below are maxima over the initial state and all four endpoints.

| Clock and profile | Velocity-only pressure error (Pa) | Paired pressure error (Pa) | Paired maximum normalized budget residual |
| --- | ---: | ---: | ---: |
| Recursive, smooth isothermal | 7.27596e-11 | 7.27596e-11 | 8.98296e-15 |
| Common, smooth isothermal | 4.36557e-11 | 4.36557e-11 | 3.24753e-14 |
| Recursive, imported nonisothermal | 6.23895e-2 | 1.27125689e-3 | 1.03161e-14 |
| Common, imported nonisothermal | 6.22812e-2 | 1.27123833e-3 | 1.24502e-14 |

The nonisothermal pressure error decreases about 49-fold with the paired
representation. The residual is almost the same for common and recursive
clocks. This supports separating the remaining thermal-contact error from
recursive clock integration; it does not prove a complete causal attribution
or establish spatial convergence. The isothermal native control already
holds pressure to roundoff with consistent velocity alone.

`native/record.json` stores exact commands, input and executable hashes,
return codes and endpoints. Per-component inventory CSVs retain physical
flux changes, numerical-energy corrections and raw residuals. The documented
budget is conservation **including** returned numerical-energy corrections,
not a claim that such corrections vanish.

The prior common-clock selected-pressure application outputs were also
compared byte for byte: all 54 CSVs match the prior integration run. See
`common-policy-parity.json` for both directories. Source hashes were checked
against `source-identity.json` after the native runs with no drift.

## Completed recursive checks

GNU Fortran 13.3.0, Debug, OpenMPI 4.1.6; the launcher is explicitly
`/usr/bin/mpiexec.openmpi` from the CMake cache. The shell's unqualified
`mpiexec` names another MPI installation and is not the tested launcher.
The existing selected chemistry configuration uses its pinned SUNDIALS build.

The focused CTest batch passed **8/8**, with zero failures or skips:

- Temporal pressure interpolation and transactional refusal unit test.
- Owner-local pressure halo and direct flux tests at 2 and 3 MPI ranks.
- Common selected-mechanism pressure application restart.
- Recursive fixed-mechanism pressure application lifecycle.
- Recursive selected-mechanism reacting/transport/regrid restart.
- Native recursive pressure saved-stage replay at 3 ranks.
- Native recursive pressure forced-refusal/retry at 3 ranks.

The initial build failed because the new `pressure_level` dummy was omitted
from a wrapped MPI halo argument list. This was corrected; `build-retry.log`
and `probe-build.log` are the successful builds used for these checks.
No test tolerances or golden references were changed to obtain a pass.

The additional disabled-pressure recursive selected application test also
passed: **9/9 distinct tests total**, zero failures or skips. Its 54 CSVs
are byte-identical to the previous velocity-only integration run. Together
with common-pressure parity this is 108 unchanged CSVs across two policies.
Provenance, project contract, Python syntax and `git diff --check` passed.
The machine-readable final record is `record.json`, including source,
configuration, binary and report hashes, toolchain, test commands and results.

## Equal-endpoint timestep refinement

Eight 25 ns steps were compared at the same 200 ns endpoint using the same
frozen binary and imported nonisothermal input. Both controls completed and
passed all component budgets (`native_dt25/record.json`).

| Clock | Paired error with 50 ns (Pa) | Paired error with 25 ns (Pa) | 25 ns maximum normalized budget residual |
| --- | ---: | ---: | ---: |
| Recursive | 1.271256886e-3 | 1.271194298e-3 | 8.25289e-15 |
| Common | 1.271238332e-3 | 1.271166155e-3 | 1.82140e-14 |

The residual changes by only 0.0049% (recursive) or 0.0057% (common) when
the timestep is halved. This is evidence that reducing the timestep alone
will not remove this residual at the tested spatial resolution. Two timestep
sizes do not establish an asymptotic order. Spatial convergence and pinned
reference qualification remain open; neither is replaced by these lifecycle
checks. No running calculation remains for this evidence batch.
