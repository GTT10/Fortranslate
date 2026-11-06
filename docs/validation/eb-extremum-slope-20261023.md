# Isolated extremum-preserving slope experiment

Source: [Sekora and Colella, arXiv:0903.4200v2, section 2.1](https://arxiv.org/html/0903.4200v2).
The scalar algorithm compares neighboring second differences near extrema;
elsewhere it applies conventional van Leer/MC limiting. The published
C_VL=1.25 is used without fitting. This is a separate linear reconstruction,
not a new acceptance test for the existing quartic candidate and not a
PeleC/AMReX default-policy implementation. The paper does not establish
positivity of an arbitrary coupled EB state.

Implemented `tools/reference_probes/eb_extremum_slope_mod.F90` and its stream
probe. Half-parent child averages are q0 minus/plus slope/4, retaining their
mean. The scalar derivation assumes regular equal-volume cells; EB weights,
multicomponent coupling, EOS admission and MPI are not connected.

At HEAD `0a0e00ebe9138145fbb0e0e745665d0d768df009` with existing dirty work,
GNU Fortran 13.3.0 Debug bounds/FPE checks completed 284 cases: 21 constant,
linear and quadratic-average controls (three scales and five extremum
positions), 41 translated steps, 102 tanh profiles and 120 sine parents.
Polynomial child averages match analytic integrals within the stated
floating-point tolerance, all child pairs preserve their parent mean,
and sampled step/tanh children remain within neighboring values. Sine
errors are measurements, not pass/fail convergence gates.

| Sine grid | EP maximum child error | Conventional MC error |
|---|---:|---:|
| 8 | 0.17908405 | 0.07417904 |
| 16 | 0.01497511 | 0.01909149 |
| 32 | 0.0003519973 | 0.004807542 |
| 64 | 0.00004426584 | 0.001204060 |

The coarse-grid regression is explicit: this is not uniformly more accurate.
These are unit-amplitude sine cell averages, not flow velocities or evolved
momentum. The fine-grid improvement justifies a saved-stencil comparison,
not production adoption or a large solver run. Next compare actual runtime
energy/momentum lines, then determine whether a conservative family with
existing EOS safeguards can be formed without reproducing the binary
admission switches. C2-C5 remain open.

Reproduce:

```sh
python3 tools/check_eb_extremum_slope.py \
  --work-dir build/eb-extremum-slope-20261023
```

The artifact directory contains inputs, outputs, compiler command and source/
binary hashes. Python syntax, mechanism provenance, project contract and
`git diff --check` pass. No time-evolution job was launched.
