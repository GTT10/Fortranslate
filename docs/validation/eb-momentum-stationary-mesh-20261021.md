# Stationary matched-candidate spatial diagnostic, 2026-10-21

The bounded-momentum candidate pair passes unchanged field/inventory gates
on 8/16/32 grids through 100 ns. Only grids8/16 were newly run (58.25/326.28 s);
the identified passing32 result was reused. Conservative analysis took12.46 s.
Both binaries use the matched private candidate, held-native energy policy,
activation-basis mechanism,129-point initial profile and50 ns steps. This
is not unmodified-PeleC parity or a rerun of the current production app.

Maximum temperature differences between the two solvers at100 ns are
1.0961e-5,7.8479e-6 and1.1258e-5 K for8,16,32. However, spatial convergence
of time evolution is not established. The apparent approximately1.997 order
of endpoint density/energy/temperature is dominated by initial restriction
error. For example, temperature RMS pair differences change from0.658840/
0.165106 K initially to0.658777/0.165086 K at100 ns.

| Evolved difference | RMS8-to-16 | RMS16-to-32 | Observed order |
| --- | ---: | ---: | ---: |
| Density (kg/m3) | 6.6612e-9 | 9.9215e-9 | -0.575 |
| Total energy (J/m3) | 0.090207 | 0.065649 | 0.458 |
| Temperature (K) | 4.5388e-4 | 3.7804e-4 | 0.264 |
| x-momentum (kg/m2/s) | 3.5308e-6 | 1.9860e-6 | 0.830 |
| y-momentum (kg/m2/s) | 1.4162e-6 | 1.5840e-6 | -0.162 |
| z-momentum (kg/m2/s) | 1.9485e-6 | 1.9450e-6 | 0.003 |

These are fluid-volume RMS values of the fine-to-coarse restricted change
minus the coarse change. Temperature uses the stored fluid-volume mean,
not EOS reconstruction. Native-held orders closely reproduce these values.
No observed-order threshold was weakened or treated as a passed gate.

A disjoint geometric partition of the16-to-32 evolved squared error assigns
74.89% of energy and75.73% of temperature to a two-cell coarse-fine band,
24.77%/24.17% to cut cells (which take priority), and under0.1% to the remaining
interior. Periodic x/y fine lower faces are coarse-fine interfaces too.
These squared-error shares are localization, not causal attribution.

An initial ideal-gas pressure audit using stored temperature and species
finds valid coarse/fine pressure uniform at101325 Pa to roundoff. Covered
coarse cells reach101327.0965,101325.5910,101325.1517 Pa on8,16,32 grids.
The conserved restriction/EOS interaction is a candidate explanation for
interface disturbances; changing pressure without conserving energy is not
an authorized numerical fix. A candidate-specific8-grid time-refinement
pilot is the next discriminator before broader calculations.

Exact hashes, commands, reports and unchanged reused32 identity are in
`build/eb-momentum-stationary-mesh-20261021/record.json`; norms in
`mesh-comparison.json`, geometric bins/script in `error-localization.json`
and `localize.py`, and pressure values in `initial-pressure-audit.json`.

## Candidate-specific time pilot

The fixed8-grid50/25/12.5 ns runs through100 ns all pass the matched field
gates. Observed temporal orders remain approximately1 in both solvers
(Fortran density1.043,energy1.036,temperature1.026,x-momentum1.009).
The50-to-25 ns RMS differences are4.405e-6 K in temperature,
0.0013903 J/m3 in energy and8.632e-7 kg/m2/s in x-momentum.
Compared with the evolved8/16 spatial differences, time-step sensitivity
is smaller for temperature/energy but material for x-momentum. This8-grid
pilot cannot establish time convergence at16/32. One16-grid25 ns run is
planned to compare the spatial8/16 pair directly at the smaller step.
Exact run identities and norms: `build/eb-momentum-time-pilot-20261021/`.

## Direct half-step spatial comparison

The new16-grid25 ns run and reused8-grid25 ns run pass all field gates.
At100 ns, halving the time step changes the evolved8/16 RMS spatial difference
by -0.390% for temperature, -0.983% for energy, +4.349% for density,
+4.742% for x-momentum, +1.351% for y-momentum and +1.520% for z-momentum.
Native-held changes agree closely. Thus halving the step does not remove
these spatial differences; it does not establish asymptotic convergence or
prove a unique interface cause. No32-grid time rerun is justified by this
result alone. Preserve the original field gates and investigate interface/
initialization consistency before changing numerical policy further.
Exact run identities, norms and ratios:
`build/eb-momentum-halfstep-mesh-20261021/{record,mesh-comparison,step_sensitivity}.json`.
