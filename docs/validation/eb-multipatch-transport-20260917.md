# Multiple-patch EB molecular transport — 2026-09-17

Implementation: `ee6343fcedfcc484b3ccdc52e14e85e0ef7feb57`, based on
`0bbce08` on local `main`.

`advance_reactive_eb_multipatch_transport_3d` connects coarse and multiple fine
Euler molecular-transport updates through the existing multiple-patch lifecycle.
Coarse/fine and same-level donor halos feed each patch. Physical transport faces
remain closed; internal patch faces use donor states. Accepted limited fluxes
are reconciled by canonical box-origin ownership, restricted and refluxed,
then all refined parents are averaged. The returned limiter is the minimum of
all patch-local accepted limiters before peer reconciliation. EOS checks still
reject a nonphysical peer correction; that minimum is not a positivity proof
for reconciliation itself. State, physical-boundary inventory and limiter are
published together. The hydro wrapper keeps the same shared lifecycle.

Scope: serial, static two-level boxes, common interval, optional viscosity,
conduction, species diffusion and barodiffusion. This does not connect a
multiple-patch reacting composition, subcycling, application/restart or MPI.
Canonical flux ownership gives conservation, not general equivalence of viscous
transverse stencils to a merged patch; that broader stencil qualification remains
open.

## Independently executed checks

GNU Fortran 13.3.0, Debug/Ninja, existing
`build/selected-tree-mpi-debug`; the expanded unit is serial.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug \
  --target test_eb_multipatch_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug \
  -R '^unit_eb_multipatch_3d$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Final **1/1 passed**, zero failures/skips, 2.19 s (2.22 s total).
Provenance, contract and whitespace checks passed; no Python changed.
Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

The same unit now runs hydro and transport for refinement ratios 2 and 3,
with complete shared faces, partial faces and separated patches touching the
physical domain. Three 1e-5 s transport intervals check all-component composite
conservation, parent averaging, nonzero evolution, limiter range and zero
physical-boundary exchange. The partial-face transport case includes quadratic
transverse velocities to exercise viscous stress/work. Reversing patch order
reproduces all states, temperatures, ledgers and limiter exactly. A zero-velocity
rectangular split region agrees with a joined transport patch within the existing
`3e-12 * max(1,abs(value))` bound. This joined comparison does not establish
transverse viscous-stencil parity. Invalid overlap and a large nonphysical step
preserve all state, inventory and limiter outputs.

Logs: `build/eb-multipatch-transport-checks/{build,ctest,final-build,final-ctest}.log`.
The final log owns the reported counts; the earlier successful run preceded the
transverse-velocity fixture. No broader matrix, installed audit, external solver
parity, convergence study or physical validation was rerun.
