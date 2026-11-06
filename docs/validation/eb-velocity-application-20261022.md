# Velocity-consistent interpolation in ordinary EB-tree applications, 2026-10-22

The ordinary gas EB-tree application now accepts `consistent_velocity=true`
with `separate_momentum=true`. The selected core interpolation uses the donor
velocity times the limited density slope plus a limited relative-momentum
slope. This connects the previously diagnostic-only velocity candidate to
production application dispatch; the separate paired pressure-reconstruction
candidate remains diagnostic-only.

## Integration boundary

Serial and MPI hydro/transport prepared halos receive the new option. Common
and recursive held-MOL stages forward it through retries and replay calls.
Stored intervals and native source prefixes retain the option, and reuse with
a different selection is refused. MPI entry points compare both policy bits
collectively; the velocity option requires separate momentum. The existing
held MOL, no redistribution, no enabled LES and no spray restrictions remain.

The application checkpoint envelope gains an explicit velocity-policy marker
only when enabled. Old/default files remain unchanged when disabled. Enabling
or disabling the option across restart is refused before state/output
publication. Mixed-rank input policy differences are refused collectively.

This does not claim machine-precision pressure equilibrium in the ordinary
application: the earlier equilibrium demonstration combined this velocity
candidate with a distinct diagnostic pressure representation. Nor does it
close nonisothermal reference-field limits or physical convergence.

## Completed selected-mechanism checks

The selected H2/O2/N2 application passed common and recursive hot CVODE
restart tests (2/2, zero failures/skips). They include serial/MPI field and
snapshot comparisons, reacting and reaction-off cases, molecular transport,
topology collapse/recreation, rebalance, two-to-three-rank restart, checkpoint
policy refusal in both directions and mixed-rank policy refusal. Endpoint is
300 ns; the tests retain their existing tolerances.

| Clock | Serial maximum budget residual | Hydrogen reacted fraction |
| --- | ---: | ---: |
| Common | 7.30392e-16 | 2.29522552e-4 |
| Recursive | 5.32666e-16 | 2.29522631e-4 |

The two CTests took 53.31 s and 90.74 s (144.09 s total). This is bounded
reacting application evidence, not ignition validation or original-PeleC parity.

## Fixed-mechanism and existing-policy checks

Fixed-mechanism ordinary serial/MPI executables also pass the common and
recursive lifecycle tests (2/2; 45.72 and 63.95 s, 109.72 s total). These check
field/inventory agreement, topology changes, two-to-three-rank continuation,
policy mismatch in both directions, unsupported configuration refusal and
preservation of existing outputs on late failed restart. The legacy non-held
smoke explicitly disables both optional momentum flags.

Two existing selected-mechanism tests also pass: converged-EOS default-policy
CVODE restart and the original separate-momentum common restart (45.66 and
56.49 s, 102.19 s total). Across these three runs, six CTests passed with zero
failures or skips. These counts do not include earlier build failures.

The ordinary namelist selection is:

```fortran
separate_momentum = .true.
consistent_velocity = .true.
```

Use it with `time_integrator='sdc'`, `mol_mode=true`, `dual_energy_enabled=true`
and `redistribution='none'`, without enabled LES or spray. The defaults remain
false. The fixed binaries tested are `pelef_reactive_eb_tree_3d` and
`pelef_mpi_reactive_eb_tree_3d`; selected-mechanism tests compile the same
ordinary application source with the selected H2/O2/N2 runtime.

## Provenance

Base HEAD `0a0e00ebe9138145fbb0e0e745665d0d768df009`, dirty main. GNU Fortran
13.3.0 Debug, OpenMPI in the existing isolated registration build, configured
SUNDIALS and pinned h2o2_full selected mechanism. No compiler/MPI change.
Artifacts: `build/eb-velocity-application-20261022/`, with build/test logs,
source identities and linked application reports. Earlier build failures
(extra MERGE argument during forwarding and a stale preprocessed multiline
condition) were corrected; failed logs remain and are not counted as passes.
Existing comparison/conversion and executable-stack warnings remain.

The application test scripts now support `--consistent-velocity`; the selected
and fixed common/recursive variants are registered CTests. Provenance, Python
syntax, project contract and whitespace checks pass. A full configuration
matrix and installed-binary audit are not claimed.

Exact completed build targets and test selections:

```sh
cmake --build build/eb-eos-application-20261021/isolated-openmpi/registration --target test_selected_full_reactive_eb_tree_3d test_selected_full_mpi_reactive_eb_tree_3d --parallel 2
cmake --build build/eb-eos-application-20261021/isolated-openmpi/registration --target pelef_reactive_eb_tree_3d pelef_mpi_reactive_eb_tree_3d --parallel 2
ctest --test-dir build/eb-eos-application-20261021/isolated-openmpi/registration -R '^regression_selected_eb_tree_velocity_restart_' --output-on-failure -j 1
ctest --test-dir build/eb-eos-application-20261021/isolated-openmpi/registration -R '^regression_eb_tree_velocity_application_' --output-on-failure -j 1
ctest --test-dir build/eb-eos-application-20261021/isolated-openmpi/registration -R '^regression_selected_eb_tree_(converged_cvode_restart|momentum_restart_common)$' --output-on-failure -j 1
```

`record.json` lists the four tested executable hashes and application-report
hashes. Source identity is recorded separately in `source-identity.json`.
