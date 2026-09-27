# Short decane jet connection

`tools/prepare_jet_reference.py` connects the pinned Spray-Jet material and
rate data to the selected MPI EB tree spray application. It checks the original
input and ROI SHA-256 values plus decane YAML/Simple-fit provenance before
creating a new output directory. Original reference files are preserved there.
The SI rate file contains all 7,502 upstream knots. External rate files in the
MPI EB tree application accept up to 16,384 points; inline namelist profiles
retain their 64-point limit. Restart binds every numerical knot, including
those outside the current run interval. The event schedule still controls
which part of the table is injected; this case uses two active intervals.

```sh
python3 tools/prepare_jet_reference.py --bundle build/decane-3sp-simple-reference/decane_3sp.json --reference-source build/pelec-pinned-hit-source --output build/decane-jet-connection
```

First prepare the decane bundle as described in the
[import record](../../docs/validation/decane-spray-import-20260927.md).
Use a separate selected-mechanism build:

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/decane-simple-openmpi-debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DCMAKE_BUILD_TYPE=Debug -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_SPRAY=ON -DPELEF_ENABLE_TESTS=OFF -DPELEF_ENABLE_SUNDIALS=OFF -DPELEF_MECHANISM_BUNDLE="$PWD/build/decane-3sp-simple-reference/decane_3sp.json" -DPELEF_MECHANISM_SOURCE="$PWD/build/decane-3sp-simple-reference/decane_3sp.yaml"
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/decane-simple-openmpi-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
# From the prepared output directory:
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none /usr/bin/mpiexec -n 1 ../decane-simple-openmpi-debug/pelef_mpi_reactive_eb_tree_3d_selected_spray tree.nml spray.nml
```

The source uses nitrogen at 440 K and 6 MPa, 363 K decane, 4.5 micrometre drops,
a 90 micrometre nozzle and discharge coefficient 0.89. The 21 degree full
spread angle becomes a 10.5 degree half-angle. ROI times and rates convert
from ms/g/s to s/kg/s. The first positive knot starts injection at 0.2
microseconds; two event intervals end at 0.6 microseconds. Their integrated
masses and average-rate velocities are recorded in `preparation.json`.

This is an explicit short adaptation: 4x8x4 all-fluid cells, no AMR, four
deterministic parcels per event, and a nozzle center displaced 1e-8 m inside
the lower-y wall. Fuel reference temperature is 300 K because upstream's
298.15 K lies below the NC10H22 thermo range; the source thermo bounds are
unchanged. The Fortran case uses interval-integrated mass and event-average
velocity. It does not reproduce upstream's stochastic injection trajectory,
full 3.8 ms jet evolution or refined mesh. Chemistry is OFF.
