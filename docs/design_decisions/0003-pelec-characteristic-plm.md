# Decision 0003: separate time-traced PeleC PLM from the SSPRK2 PLM path

## Status

Accepted for the one-dimensional constant-`gamma` milestone.

## Context

The existing `plm` mode reconstructs componentwise primitive face values and uses SSPRK2 to obtain second-order time accuracy. PeleC `Source/PLM.H` instead projects slopes into characteristic waves and traces those waves over `dt/dx`, producing time-centered interface states before the Riemann solve.

Applying SSPRK2 unchanged to those already time-centered face states would mix two distinct time discretizations and obscure parity with the source algorithm.

## Decision

Add a separate `pelec_plm` reconstruction mode that:

- retains the tested primitive MC/minmod slope calculation;
- projects the slope into left acoustic, contact, two shear, and right acoustic waves;
- performs the one-dimensional `u-c`, `u`, `u+c` tracing used by the relevant PeleC PLM kernel;
- sends time-centered face states to either selectable Riemann solver;
- advances with one conservative Godunov update;
- leaves `pcm` and componentwise `plm` on their existing SSPRK2 path.

## Consequences

Benefits:

- reconstruction and flux changes remain independently selectable;
- the original robust baselines remain available;
- characteristic algebra and time tracing have dedicated unit tests;
- conservation is maintained by one shared face flux per interface.

Limitations:

- slopes are currently second-order MC/minmod, not PeleC fourth-order slopes;
- flattening is absent;
- the thermodynamic state is constant-`gamma` ideal gas;
- general-EOS internal-energy, species/passive scalars, EB, and multidimensional corrections are absent.

The next step is fourth-order limited slopes plus flattening, followed by a strong-shock Sedov gate.
