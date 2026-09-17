# Plain Euler EB patch diffusion timestep — 2026-09-17

Implementation: `71316253af8d53b7cc46dc817d5b3130d387e3c0`, based on
`bd2b39f` on local `main`.

`reactive_eb_transport_patch_timestep_3d` reuses the existing thermodynamic and
transport evaluation to obtain the largest enabled momentum, thermal or
species diffusivity. It then estimates an Euler patch interval from
`transport_cfl / (Dmax * max_cell(metric))`, with
`metric = sum_d((alpha_minus + alpha_plus)/(2*h_d**2))/kappa`.
The allowed transport CFL remains `(0, 0.5]`.

Optional `donor_diffusivity` raises the diffusivity maximum for remote face
donors. Callers must supply that bound when appropriate. Outer open faces
remain in the estimate even when physical transport flux is disabled, making
that part conservative. Embedded walls have zero molecular flux and no added
conductance. All-disabled transport returns an unlimited interval and zero
diffusivity. Invalid inputs reset both scalar outputs. Existing redistributed
transport timestep arithmetic is unchanged.

This is a geometry-scaled diffusion estimate, not a positivity or stability
proof for variable-property multicomponent transport, barodiffusion, reflux or
reaction coupling. It does not replace the accepted-flux species limiter or
transactional rejection. Application/hierarchy timestep selection is still
not connected by this increment.

## Independently executed checks

GNU Fortran 13.3.0, Debug/Ninja, existing
`build/selected-tree-mpi-debug`. The exercised unit is serial.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug \
  --target test_reactive_eb_transport_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug \
  -R '^unit_reactive_eb_transport_3d$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Final **1/1 passed**, zero failures/skips, 0.91 s (0.94 s total).
Provenance, project contract and whitespace checks passed; no Python changed.
Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

The unit preserves existing transport checks and adds an anisotropic regular
patch comparison against the existing interval. In each Cartesian orientation,
a planar cell with fluid fraction 0.05 has metric `12/h**2`, so its interval is
one quarter of the regular-grid estimate. Doubling remote diffusivity halves
that interval. Disabled transport, ignored covered storage, negative donor
bounds and invalid active states exercise the corresponding output contracts.

Logs: `build/eb-patch-diffusion-cfl-checks/{build,ctest,final-build,final-ctest}.log`.
The final log owns the reported counts; the earlier successful run preceded the
regular-patch check and overflow guards. No broader matrix, installed audit,
hierarchy stepping, external solver parity or physical validation was rerun.
