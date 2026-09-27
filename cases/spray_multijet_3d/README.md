# Four-port decane connection

Prepare the pinned Spray-Multijet adaptation using the existing
[decane bundle](../../docs/validation/decane-spray-import-20260927.md):

```sh
python3 tools/prepare_multijet_reference.py --bundle build/decane-3sp-simple-reference/decane_3sp.json --reference-source build/pelec-pinned-hit-source --output build/decane-multijet-connection
```

Run the selected MPI EB tree spray executable from the output directory with
`tree.nml spray.nml`. Build instructions are in the
[single-jet guide](../spray_jet_3d/README.md). After completion:

```sh
python3 tools/check_jet_inventory.py --run-dir build/decane-multijet-connection --output build/multijet-inventory.json
```

Four ports occupy the x/z quarter points of a 12 mm cube, with periodic x/z,
a lower-y no-slip wall and upper-y outflow. Each injects at 600 m/s and
0.0023 kg/s. Gas starts at 900 K and 2.93 MPa with O2/N2 mass fractions
0.233/0.767. Eight scheduled events cover two 100 ns intervals; four
deterministic parcels per event give 32 parcels and 1.84e-9 kg injected fuel.

The original input, initialization code and build configuration are byte-pinned
in `references/spray_multijet_case.json` and copied into the prepared directory.
The upstream species is NC10H22, but its liquid-property comment says dodecane;
this adaptation preserves the input values: density 669 kg/m³, heat capacity
2593.9 J/(kg K), latent heat 302940 J/kg and its Antoine coefficients. It does
not substitute the different single-jet properties or certify these as decane
material data.

Explicit adaptations are a 4x8x4 all-fluid mesh without AMR, deterministic
sampling, a 1e-8 m inward nozzle displacement, 300 K fuel reference temperature
instead of the out-of-range 298.15 K, and an Antoine evaluation interval of
300–489 K. That interval is an application restriction, not a fitted validity
claim. Chemistry is OFF. Full-duration stochastic jet/reference parity,
refinement and physical validation remain open.

See the [bounded execution record](../../docs/validation/multijet-connection-20260927.md).
