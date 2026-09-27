# Pinned Spray-Evaporation reference preparation

This prepares a **3D periodic adaptation**, with 400 K gas and 310 K liquid,
of pinned PeleC's Spray-Evaporation. It does not reproduce the original
2D 273 K/272 K case with legacy `Hard` boundaries. The original temperatures
are outside the imported NC7H16 NASA7 range, which remains unchanged.

Use the pinned source export described in
[the HIT guide](../reactive_hit_3d/README.md). Particle support changes the
AMReX configuration, so use a separate build directory:

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH GIT_CEILING_DIRECTORIES="$PWD/build" cmake -S build/pelec-pinned-hit-source -B build/pelec-pinned-spray -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DPELE_ENABLE_MPI=OFF -DPELE_ENABLE_PARTICLES=ON -DBUILD_TESTING=OFF
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/pelec-pinned-spray --target PeleC-Spray-Evaporation --parallel 2
python3 tools/prepare_evaporation_reference.py --bundle build/heptane-3sp-reference/heptane_3sp.json --output build/spray-evaporation-reference-run
```

The preparation tool verifies the pinned YAML SHA and requires a new output
directory. It converts the 79/21 N2/O2 mole mixture to mass fractions using the
bundle molecular weights, and computes parcel multiplicity from 1e-10 kg total
liquid, four parcels, 20 µm diameter and 681.4 kg/m³ density. It writes a CGS
gas input, particle ASCII input and preparation metadata. Parcel fields follow
the pinned `SprayComps` ordering, with x varying first in position/velocity.

After a successful build, run `PeleC-Spray-Evaporation reference.inp` from the
prepared directory. The executable is under
`build/pelec-pinned-spray/Exec/RegTests/Spray-Evaporation/`.
The input requests two `1e-7` s steps, fourth-order PLM and two SDC iterations. Reference
gas and particle outputs will be written for initialization and each step.

## Current evidence boundary

Preparation implementation: `b9b4f9c0`. Source identities are the unchanged
four-commit export in `references/pelec_baseline.json`. Configuration completed
with GNU C/C++ 13.3.0 Release, serial CPU, particles ON. At preparation recording,
the reference build was still running; neither executable completion nor a
reference trajectory comparison is claimed. Logs are
`build/pelec-pinned-spray-configure.log` and `build/pelec-pinned-spray-build.log`.
Input preparation, Python syntax, project contract and diff checks passed.

The earlier [Fortran integration fixture](../spray_3d/heptane_nonreacting.nml)
has matching basic material/initial scales but uses the serial spray
application's different gas/particle integration path. Scheme, staging,
deposition and property conventions must be aligned before treating a field
comparison as solver parity. No tolerance or reference material data has been
changed to hide those differences.

Preparation now also creates `tree.nml`, `spray.nml` and an all-fluid nodal
`level_set.dat` for the selected MPI EB tree particle SDC application. This
uses the source molecular weights to set gas density at 400 K. Build it in a
separate MPI-enabled selected-mechanism configuration and run
`mpiexec -n 1 pelef_mpi_reactive_eb_tree_3d_selected_spray tree.nml spray.nml`.
The one-rank integration check is recorded
[here](../../docs/validation/nonreacting-particle-sdc-20260927.md); full reference
trajectory parity remains pending.
