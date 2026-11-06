# Paired pressure in ordinary EB-tree applications, 2026-10-22

`paired_pressure=true` is connected to ordinary serial/MPI common-clock gas
EB-tree applications. It requires both velocity-consistent and separate
momentum interpolation, held MOL without redistribution, and no rotation.
The defaults retain the previous policy. Recursive subcycling is explicitly
refused; it requires a distinct time-consistent pressure-history construction.
See the [numerical decision](../design_decisions/eb-paired-pressure.md).

The implementation builds EOS pressure before fluid-volume restriction,
exchanges the restricted pressure stencil, and reconstructs both fine halos
and covered coarse states. It retains conserved inputs and uses transactional
EOS admission. Transport and hydro share this representation. Saved stages,
replay policy checks and checkpoint envelopes include the selection.

Artifacts, exact commands, source/binary hashes and completed results are in
`build/eb-pressure-application-20261022/`. This is application integration and
focused conservation/rollback evidence, not a full configuration matrix,
installed-binary audit, original-PeleC parity or physical validation.

## Completed checks

Seven distinct CTests have a completed passing result (zero final failures
or skips): transfer kernel; two/three-rank halo tests; fixed common-clock
pressure application; selected hot-CVODE pressure restart; existing selected
converged-EOS default and velocity-recursive application checks.

The new halo test uses isobaric nonuniform-temperature N2 on three tree levels,
two ownership layouts and two/three ranks. It checks serial/MPI gas, held
energy, temperature and volume agreement, unchanged input gas/held fields,
pressure within 1e-7 Pa of 60000 Pa, interface velocity within 1e-9 m/s of
10 m/s, and the Riemann momentum-flux identity `F_mx - 10 F_rho = 60000`
to 1e-7 Pa. Mixed policy and remote invalid held energy leave every output
halo unchanged. These are initial spatial/flux checks, not a time-evolution
accuracy claim.

Fixed and selected application tests cover chemistry/transport, topology
collapse and recreation, owner migration, 2-to-3-rank restart, policy mismatch
in both directions, and preservation of checkpoint/output data on refusal.
The selected 300 ns reacting case has maximum serial budget residual
9.22354e-16, MPI/restarted residual 2.37695e-16, and hydrogen reacted fraction
2.29522552e-4. The fixed application also checks refusal of recursive pressure,
missing velocity consistency, rotation, missing held energy, non-MOL stepping
and redistribution. No full forced-retry/replay suite was rerun for the new
pressure policy; stored identity is wired and checkpoint/halo rollback is
covered by these checks.

Existing-policy output parity is exact: all 46 default-policy CSV artifacts
and 54 velocity-recursive CSV artifacts match the previous independently
completed application executions byte for byte. Mechanism provenance, Python
syntax, project contract and whitespace checks pass. No tolerances or golden
hashes changed.

## Failures retained and corrected

The first build failed because the mechanical argument threading appended a
pressure argument inside three `MERGE` calls. Corrected calls and independent
pressure-policy bits compile successfully. The first two halo tests failed
while creating their new fixture: its inherited light mixture made the chosen
density/pressure fall below the EOS temperature range. The fixture now uses
explicit N2 (approximately 390--620 K); both rank counts pass, including the
added direct Riemann and rollback assertions. Original failure logs are
retained separately from the passing rerun. The remaining five CTests passed
on their first execution (0.07, 91.41, 44.87, 54.81 and 59.09 s); corrected
halo tests passed in 8.43 and 6.99 s.

Exact source identity is the dirty `main` tree based on
`0a0e00ebe9138145fbb0e0e745665d0d768df009`, with full source SHA256 inventory
and executable hashes in `record.json` / `source-identity.json`. Build and
test argv, GNU/OpenMPI versions, CMake configuration hash, application report
hashes and all logs are retained in the artifact directory.
