# Evaporation temperature inversion audit — 2026-09-27

Audited Fortran source: `595bc61e13dde5659429a1976a8a23738ef8a623`.
Pinned PeleC: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`;
PelePhysics: `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`.
This audit changes no numerical implementation or acceptance tolerance.

## Finding

The 64-cell, two-step evaporation reference stores temperatures differing
from temperatures reconstructed from conserved energy by up to
`8.872212333699281e-7 K`. Cantera 3.2.0, using the pinned heptane YAML,
reconstructs the corresponding Fortran temperatures within `1.514e-10 K`.
The reference discrepancy is the scale of the remaining stored-temperature
comparison residual.

The pinned `PeleCUtilities.H::pc_cmpTemp` uses the stored temperature as the
initial guess for `Fuego::REY2T`, which calls `GET_T_GIVEN_EY`. In the pinned
heptane mechanism, the default iteration stops **before applying** a Newton
correction smaller than `1e-6 K`. With the default small/negative-energy
settings, `pc_rst_int_e` obtains internal energy from total energy minus
kinetic energy; this case does not require a separate dual-energy explanation.

## Focused compiled check

Using GNU C++ 13.3.0 and the existing pinned reference include/define flags,
a standalone harness calls the actual mechanism function on all 64 reference
CSV states. Specific internal energy and mass fractions are reconstructed
directly in CGS from those states. Results:

| Kernel configuration | Maximum change from stored T | Nonzero error codes |
| --- | ---: | ---: |
| Default | 0 K | 0 |
| `CONVERGENCE` (`1e-12 K`) | `8.872212333699281e-7 K` | 0 |

The maximum initial Newton correction is `8.872212394557612e-7 K`.
Thus the compiled default inversion retains every stored temperature, whereas
the tighter diagnostic kernel recovers the independently reconstructed value.
Only the standalone diagnostic uses `CONVERGENCE`; the reference application
and its output are unchanged.

Artifacts: `build/spray-temperature-kernel-audit/` contains the harness,
64 input states, both compiled executables and outputs, exact commands, source
identity and CSV SHA-256 in `comparison.json`. The independent Cantera audit
is `build/spray-face-transport/temperature-eos-audit.json`. Input reference CSV:
`build/spray-evaporation-euler/reference-gas.csv`.

This explains why stored temperature is not an exact common EOS observable
for this comparison. It does not assign the entire inter-solver energy or
momentum difference to that stopping rule: no full application intervention
with a tighter reference inversion was run. Preserve the Fortran conservative
energy update. Full spray parity and physical validation remain open.

Mechanism provenance, Python syntax for the contract checker, project contract
and `git diff --check` pass. No solver rebuild, full matrix, restart rerun or
physical validation was needed for this source/kernel audit.
