# Selected single-level 3D EB

`pelef_reactive_eb_3d_selected` uses the configured mechanism bundle for planar,
extruded and nodal level-set geometry. Its shared application now integrates
source element composition, including carbon and inert elements, rather than
recognizing a fixed list of species labels. The reaction integrator remains
the bundle's configured integrator.

`methanol.nml` is a small uniform-density 21-species FFCM1_Red case for the
bundle prepared by `tools/prepare_spray_reference.py`; use that imported JSON
as `PELEF_MECHANISM_BUNDLE`. `generic_element.nml` exercises the synthetic A/B,
X-element test bundle; it is a software fixture, not a physical mechanism.

Selected checkpoints use schema 3, binding source element names and atom
counts along with the bundle SHA, integrator and composition. Earlier schema-2
files are not migrated by this frontend. Fixed-runtime callers keep schema 1;
library callers omitting generic metadata retain their legacy three-element
path. Checkpoint initial element amounts are now populated even for selected
runs with chemistry disabled.

The shared application still gates global invariant and elemental inventories
without a general external-boundary exchange ledger. The qualified methanol
fixture has uniform density; a density-contrast trial exceeded the unchanged
mass gate and remains unqualified. This is not an ignition or transport-accuracy
validation. See [evidence and limitations](../../docs/validation/eb-single-level-elements-20260917.md).
