# Selected single-level 3D EB

`pelef_reactive_eb_3d_selected` uses the configured mechanism bundle for planar,
extruded and nodal level-set geometry. Its shared application now integrates
source element composition, including carbon and inert elements, rather than
recognizing a fixed list of species labels. The reaction integrator remains
the bundle's configured integrator.

`methanol.nml` is a small density-contrast 21-species FFCM1_Red case for the
bundle prepared by `tools/prepare_spray_reference.py`; use that imported JSON
as `PELEF_MECHANISM_BUNDLE`. `generic_element.nml` exercises the synthetic A/B,
X-element test bundle; it is a software fixture, not a physical mechanism.

Selected checkpoints use schema 6, binding source element names and atom
counts along with the bundle SHA, integrator and composition. Earlier schema-2/3
files are not migrated by this frontend. Fixed-runtime applications use schema 4;
library callers omitting exchange vectors retain schemas 1/2/3. Checkpoint initial element amounts are now populated even for selected
runs with chemistry disabled.

The shared application accumulates actual hydro boundary and wall exchanges
and chemical species changes, retains them through restart, and checks all
component budgets. Element checks subtract external exchange only. Diagnostics
print initial/final amounts, each exchange and the raw residual. Restart
budgets are checked before continuation; checkpoint output is checked before
writing. The original relative gates are retained, with a step-scaled floating
point roundoff allowance for species components that were formerly excluded
when chemistry was enabled.

The density-contrast methanol fixture now passes the boundary-corrected mass
budget and exact state/budget continuation. This is not an ignition or
transport-accuracy validation. See [evidence and limitations](../../docs/validation/eb-application-exchange-20260917.md).
