# Pinned Spray-EB input connection

Generate the upstream 64³, 10 cm cube, 1500 K gas/wall case and its five
initial heptane droplets:

```sh
python3 tools/prepare_spray_eb_reference.py \
  --bundle build/heptane-3sp-simple-reference/heptane_3sp.json \
  --reference-source build/pelec-pinned-hit-source \
  --output build/spray-eb-reference-inputs
```

The source export must match `references/spray_eb_case.json`. The bundle must
be the pinned reaction-free heptane mechanism with regenerated Simple fits.
The output includes original source copies, SI namelists, the sign-converted
nodal plane and `preparation.json`. Existing output directories are rejected.
Use the selected MPI EB spray executable built against that bundle. Execute
from the generated directory so the relative level-set path resolves.

The default preserves the original 64³ grid and final time
`0.000273365690540257 s`. `--cells` and `--final-time` can select a shorter
connection check; every change is recorded. No execution is launched by the
preparer. The five records preserve position, velocity, temperature, diameter,
composition and one-real-droplet multiplicity, represented as five t=0 events.
The ASCII schema comes from pinned PelePhysics `SprayFuelData.H` and
`SpraySetup.cpp`: xyz then velocity, temperature, diameter, composition,
number density, initial count, two breakup variables, and film height.

The generated input now selects `mol_mode=.true.` in the coupled integrator,
matching the upstream choice of method-of-lines integration. The Fortran MOL
reconstruction adds a positive-face slope safeguard beyond the pinned formulas.
Other unresolved differences include a single root patch and a 100 ns upper
timestep. The Antoine evaluation range is restricted to
300–371.6 K. Particle contacts use the current specular path; full wall-impact
reference fidelity is not established. This is not full Spray-EB parity.

The 8³, 200 ns connection now completes with five particles after repairing
physical-wall source donor positions. The earlier failed particle and gas-only
artifacts remain preserved. See [repair and evidence](../../docs/validation/spray-eb-source-halo-20260927.md).
The full 64³ duration and reference-solver comparison remain unverified.

The updated 8³, 200 ns MOL connection also completes on two MPI ranks with
five particles, embedded wall heat and no-slip traction. See
[MOL connection evidence](../../docs/validation/spray-eb-mol-connection-20260927.md).
