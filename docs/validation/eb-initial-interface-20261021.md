# Initial coarse/fine pressure and flux audit, 2026-10-21

The new optional `--audit-initial [--separate-momentum]` mode of
`tools/reference_probes/eb_tree_reacting_probe.F90` calls the production held
MOL spatial operator before reaction initialization or time integration.
It reconstructs EOS pressure from the returned conserved/held-energy halos
and writes pressure and coarse/fine centroid face flux diagnostics. It requires converged V4 EOS.
The ordinary two-argument trajectory path is preserved. No restriction,
interpolation, thermodynamic or flux policy has been changed.

`tools/audit_eb_initial_interface.py` runs eight small controls and pairs
fine faces with coarse faces, including periodic lower x/y interfaces.
Flux is positive-coordinate, intensive SI, with dt divided out. Comparisons
multiply each flux by its open-area fraction and sum the four fine-face
contributions with a factor of 1/4. The z lower physical wall is excluded.
It also checks uniform stationary mass/energy flux and area agreement.

## Completed observation

The saved stationary 8-cubed V4 import from
`build/eb-momentum-stationary-mesh-20261021/grid8/matched_eb_amr_tcu56ulb/initial.dat`
was reused. Both limiter policies give identical initial results:

| Quantity | Planar cut-cell control | Regular geometry ablation |
| --- | ---: | ---: |
| Maximum adjacent halo pressure offset (Pa) | 0.275222 | 0.275222 |
| Maximum interface mass flux (kg/m2/s) | 1.18231e-4 | 1.18231e-4 |
| Pressure jump / mass flux Pearson correlation | -0.998230 | -0.999907 |
| Maximum open-area weighted mass flux mismatch (kg/m2/s) | 1.95225e-5 | 1.24194e-5 |
| Maximum open-area weighted normal momentum mismatch (Pa) | 0.127029 | 0.127029 |
| Maximum open-area weighted total-energy flux mismatch (W/m2) | 64.9799 | 38.7888 |

Valid-cell EOS pressure is uniform at 101325 Pa within 4.37e-11 Pa.
Covered coarse EOS pressure differs by up to 0.140001 Pa. This uses held
energy EOS inversion, whereas the earlier 2.0965 Pa covered-coarse audit
used stored temperature; those are different observables and must not be
substituted for each other. The adjacent halo pressure offset is larger
than the covered-coarse offset, so this is not evidence of direct copying
of the donor pressure error alone.

All four uniform controls have exactly zero interface mass flux and exactly
zero coarse/fine flux mismatch in all three reported components. All eight
controls have exactly matching coarse/fine open areas. There are 48/80
coarse interface faces and 160/320 open fine interface faces in the cut/regular
controls. The regular control changes geometry only, keeping imported states,
including the earlier covered-coarse restriction. It is not an independently
initialized periodic reference run; x/y are periodic and z retains slip walls.

These observations localize an initial spatial pressure disturbance to the
nonuniform coarse/fine state construction, present before chemistry and
without cut cells. Correlation is not a unique causal decomposition between
restriction, interpolation and EOS. The next discriminator is a conservative
restriction/interpolation stencil audit, separating EOS pressure from stored
T and identifying which stencil produces the halo pressure error. Do not
flatten pressure by overwriting energy or weaken conservation to remove it.
This does not close fixed-PeleC trajectory or spatial convergence gates.

## Verification and identity

GNU Fortran 13.3.0 Debug; existing OpenMPI registration build (no toolchain
or MPI configuration changes). The probe and two affected operator test
binaries were rebuilt. CTest `^unit_eb_(tree_halo|mol_hydro)_3d$` passed 2/2,
zero failures/skips. Python periodic-face/open-area and nonfinite-refusal
checks passed 2/2. The eight-control application diagnostic completed.
Mechanism provenance, Python syntax, project contract and `git diff --check`
passed. This is focused diagnostic evidence, not the full matrix or a native
PeleC rerun.

Artifacts, commands, source hashes, compiler identity and logs:
`build/eb-initial-interface-20261021/record.json` and `final-controls/report.json`.
The base HEAD is `0a0e00ebe9138145fbb0e0e745665d0d768df009` with existing
uncommitted numerical work; HEAD alone does not identify the tested source.

Reproduce with the selected LiDryer reference probe target:

```sh
python3 tools/audit_eb_initial_interface.py \
  --probe BUILD/tests/test_eb_tree_reacting_reference \
  --input SAVED/initial.dat --output NEW_OUTPUT_DIRECTORY
```
