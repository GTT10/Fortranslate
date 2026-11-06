# Local thermal-mixing discriminator, 2026-10-22

A new `--audit-mixing` probe mode evaluates local convex combinations of
adjacent unrefined coarse fluid cells. It invokes no spatial flux, reflux,
AMR interpolation, time integration, chemistry or transport. Covered coarse
cells and pairs crossing into the refined region are excluded. The input
states and all production numerical paths remain unchanged.

For each of 304 pairs, the right state is mixed with its left neighbor using
alpha=|u_right|*dt/dx, alpha/2 and alpha/4. Here alpha=0.0004 (10 m/s,
50 ns, 1.25 mm). Gas conserved quantities and held internal energy are
combined with the same weights; the production converged EOS recovers
pressure and temperature. A second EOS evaluation holds the mixed density
and species fixed and asks what held energy would recover the linearly
weighted donor pressure. That energy gap is reported, never applied.
This is a local mixing experiment, not the actual second-order MOL update
or a decomposition of its endpoint error.

Three variants separate effects:

- 0: retain both imported donors.
- 1: reconstruct the left donor at its original pressure and temperature,
  but with the right donor's composition.
- 2: reconstruct the left donor at its original pressure and composition,
  but with the right donor's temperature.

The control reconstructions change the hypothetical pair, not the saved
input or physical trajectory. The isothermal H2/N2 contact is also run through
all three variants. There are 304 pairs x 3 variants x 3 fractions per input.

## Completed results

Maximum absolute pressure departure from weighted donor pressure (Pa):

| Imported-field variant | alpha=0.0004 | alpha/2 | alpha/4 |
| --- | ---: | ---: | ---: |
| Original T/composition differences | 0.001004149919 | 0.000502176656 | 0.000251113743 |
| Equal composition, original T difference | 0.001498036683 | 0.000749169631 | 0.000374622614 |
| Equal T, original composition difference | 5.82e-11 | 4.37e-11 | 2.91e-11 |

The original-pair maximum is at coarse cell (7,1,3) with left/right donor
EOS temperatures 1980.865828 / 2019.134172 K. Donor pressures agree with
101325 Pa to roundoff. Recovering that pressure after mixing would require
subtracting 0.00330000586 J/m3 from held energy at full alpha. With equal
composition the corresponding energy gap is -0.00492309811 J/m3.
The isothermal contact controls retain pressure within 4.37e-11 Pa.

Thus temperature-varying conservative mixing in this EOS is sufficient to
produce a pressure departure without any coarse/fine algorithm. This does
not prove that all of the previously observed 0.001247 Pa trajectory error
comes from mixing; the local experiment uses different numerical weights
from the MOL face reconstruction. It does rule out attributing every nonzero
pressure departure to an AMR synchronization defect. A forced pressure reset
would have an explicit energy cost and must not be hidden in conservative
transport accounting.

The smooth isothermal contact equilibrium gate remains 1e-8 Pa. Existing
fixed-reference field tolerances and conservation gates remain unchanged;
the imported nonisothermal trajectory remains unqualified. Actual held-MOL
integration and spatial/time convergence, with a recorded numerical-energy
budget, are still required.

## Verification and provenance

Artifacts: `build/eb-thermal-mixing-20261022/controls/` (frozen binary/input,
record.json, two CSVs/logs, checker reports, source-identity.json).
The parent directory retains the earlier two-profile pilot without the
composition/temperature controls; only the controls directory supports the
table above. No long simulation was needed; the local controls completed
in under a second in the measured orchestration call.

The checker independently verifies finite output, equal donor pressure,
weighted pressure targets, three variants with identical pair coverage,
three halved fractions, and the equal-temperature pressure-equilibrium gate.
Thermodynamic values are calculated by the production EOS in the probe;
this is not an independent EOS implementation or native-PeleC comparison.

Base HEAD: `0a0e00ebe9138145fbb0e0e745665d0d768df009`, dirty main.
GNU Fortran 13.3.0 Debug, existing isolated OpenMPI registration build.
The existing executable-stack linker warning remains.

Completed checks: probe build; 5 new thermal-mixing Python tests plus 19
existing pressure diagnostic tests (24 passed, 0 failed/skipped); mechanism
provenance, Python syntax, project contract and `git diff --check` passed.
Exact probe/checker commands and return codes are in record.json. The new
commands are `probe input.dat prefix --audit-mixing` and the isothermal
control with fourth argument `--smooth-cut-contact`, followed by
`python3 tools/check_eb_thermal_mixing.py --input prefix_mixing.csv --output report.json`.
No full matrix, MPI, installed-binary audit or physical validation is claimed.
