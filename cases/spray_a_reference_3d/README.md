# Spray-A / KHRT connection

Prepare the pinned Lu n-dodecane bundle with Simple transport using the
[import instructions](../../docs/validation/dodecane-simple-import-20260927.md).
Then prepare the application inputs:

```sh
python3 tools/prepare_spray_a_reference.py --bundle build/dodecane-simple-reference/dodecane_lu.json --reference-source build/pelec-pinned-hit-source --output build/spray-a-case
```

Build `pelef_mpi_reactive_eb_tree_3d_selected_spray` in a clean directory with
MPI/spray enabled and `PELEF_MECHANISM_BUNDLE` and `PELEF_MECHANISM_SOURCE`
pointing to that bundle and its pinned YAML. Run `tree.nml spray.nml` from
the prepared directory. The [execution record](../../docs/validation/spray-a-connection-20260927.md)
contains the exact configure and run commands. Check completed output with:

```sh
python3 tools/check_jet_inventory.py --run-dir build/spray-a-case --fuel-species NC12H26 --output build/spray-a-case/inventory.json
```

Default preparation retains the physical 12×12×60 mm domain, 900 K/6 MPa gas,
N2/H2O/CO2 mole fractions, 90 µm nozzle, 89 µm injected drops, 363 K liquid,
20-degree full cone angle, discharge coefficient 0.89, polynomial density,
inverse-temperature viscosity, Antoine pressure and KHRT constants from the
pinned `Spray-A-Wbreakup/example.inp`. The referenced liquid conductivity is
not used by this free-particle fixture; splash/film physics is disabled upstream.

Explicit adaptations: a 4×4×8 all-fluid mesh without AMR, four deterministic
parcels per event, a 10 nm inward nozzle displacement, 300 K liquid reference
temperature (the gas mechanism starts at 300 K), and property evaluation bounded
to 300–658 K. These bounds are not a fit-validity claim. Chemistry is disabled,
as in the reference example. The pressure-dependent boiling cap uses normal
boiling temperature 489 K and critical temperature 658 K.

The default first two active rate-table intervals end at 0.6 µs. `--events N`
prepares more intervals, but does not reproduce the upstream 3 ms stop time,
AMR or stochastic parcel sampling. For a run with KHRT child release, add `--allow-khrt-children` to the inventory
check. This requires KHRT preparation metadata and the application's reported
particle-ID watermark and computational split count. It checks ID bounds and
fuel inventory while allowing population growth; it does not establish
breakup-dynamics parity. Full reference comparison and experimental penetration
validation remain open.
