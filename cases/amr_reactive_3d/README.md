# Static two-level reactive 3D AMR

Run the public ratio-two, full-H2O2 entropy-wave case from a build directory:

```sh
./pelef_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_entropy_wave.nml
```

The 8-by-8-by-8 periodic coarse grid refines coarse cells 3 through 6 in all
three directions. The fine patch takes two SSPRK2 substeps per coarse step,
uses linearly time-interpolated coarse boundary states, replaces all six
coarse/fine interface fluxes through reflux, and averages the covered cells
back to the coarse level. The application writes separate coarse and fine CSV
files and rejects a result that violates composite conservation, species
closure, physicality, or average-down synchronization.

The checkpoint/restart parity sequence uses the same physical case:

```sh
./pelef_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_entropy_wave_checkpoint.nml
./pelef_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_entropy_wave_restart.nml
python3 ../tools/compare_amr_reactive_3d_restart.py \
  --reference-coarse amr_reactive_3d_coarse.csv \
  --reference-fine amr_reactive_3d_fine.csv \
  --restart-coarse amr_reactive_3d_restart_coarse.csv \
  --restart-fine amr_reactive_3d_restart_fine.csv
```

The first command stops after coarse step four and writes a versioned
checkpoint. The second restores both levels and the original conservation and
reflux history, completes the requested final time, and must reproduce both
uninterrupted CSV files byte for byte. A checkpoint is accepted only when its
schema, full NASA7 species data, solver, boundary mode, mesh, patch, CFL,
physics switches, array extents, physical state, temperature recovery, and
coarse/fine synchronization match the requested run.

The rank-local sparse application accepts the same namelists and an output
prefix. PCM rank-count and changed-rank restart parity can be reproduced with:

```sh
mpiexec -n 1 ./pelef_mpi_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_entropy_wave.nml mpi_amr_np1
mpiexec -n 8 ./pelef_mpi_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_entropy_wave.nml mpi_amr_np8
mpiexec -n 2 ./pelef_mpi_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_entropy_wave_checkpoint.nml checkpoint_np2
mpiexec -n 4 ./pelef_mpi_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_entropy_wave_restart.nml restart_np4
python3 ../tools/compare_mpi_amr_reactive_3d.py \
  --reference-prefix mpi_amr_np1 \
  --candidate-prefix mpi_amr_np8 \
  --candidate-prefix restart_np4 \
  --checkpoint amr_reactive_3d.chk
```

The characteristic-PLM path uses the corresponding `*_plm*.nml` inputs:

```sh
./pelef_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_entropy_wave_plm.nml
mpiexec -n 8 ./pelef_mpi_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_entropy_wave_plm.nml mpi_amr_plm_np8
mpiexec -n 2 ./pelef_mpi_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_entropy_wave_plm_checkpoint.nml \
  checkpoint_plm_np2
mpiexec -n 8 ./pelef_mpi_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_entropy_wave_plm_restart.nml restart_plm_np8
python3 ../tools/compare_mpi_amr_reactive_3d.py \
  --reference-prefix amr_reactive_plm_3d \
  --candidate-prefix mpi_amr_plm_np8 \
  --candidate-prefix restart_plm_np8 \
  --checkpoint amr_reactive_plm_3d.chk \
  --reconstruction characteristic_plm --limiter mc
```

Coarse and fine state and temperature are rank-local throughout stepping.
Fine ownership follows the parent coarse slab, so some ranks may own no fine
cells. Checkpoint and CSV access is selected-root-only; the checkpoint payload
contains no owner map and may therefore be resumed with a different valid
rank count. Root gather/scatter is compatibility I/O rather than scalable
distributed output.

The fixed full-H2/O2 chemistry boundary uses paired inert and reacting
hotspots:

```sh
./pelef_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_hotspot_control.nml
./pelef_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_hotspot_chemistry.nml
mpiexec -n 8 ./pelef_mpi_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_hotspot_chemistry.nml chemistry_np8
python3 ../tools/check_amr_reactive_chemistry_3d.py \
  --control-prefix amr_hotspot_control_3d \
  --reactive-prefix amr_hotspot_chemistry_3d \
  --nx 8 --ny 8 --nz 8 \
  --i-lower 3 --i-upper 6 --j-lower 3 --j-upper 6 \
  --k-lower 3 --k-upper 6 --ratio 2
python3 ../tools/compare_mpi_amr_reactive_3d.py \
  --reference-prefix amr_hotspot_chemistry_3d \
  --candidate-prefix chemistry_np8 --final-time 1.0e-9 \
  --reconstruction characteristic_plm --limiter mc
```

Each chemistry phase advances both levels privately, averages down, and
recovers the covered coarse temperature. The application gates Euler and
H/O/N conservation separately from reacting species activity. A rank with no
fine cells still participates in reaction/integrator and success consensus.
Chemistry-enabled checkpoint/restart is rejected because schema 2 does not
store its reaction and source-integrator context.

This case qualifies serial and rank-local sparse, static, strictly interior,
single-patch PCM or characteristic-PLM hydrodynamic AMR with optional fixed
elementary/full-H2/O2 chemistry and periodic coarse boundaries. Selected
mechanisms, chemistry restart, molecular transport, CTU/PPM, dynamic
regridding, scalable I/O, embedded boundaries, and boundary-touching patches
remain outside this case.
