# Sparse MPI AMR reactive hotspot

This input drives the public `pelef_mpi_amr_reactive_1d` application through
solution-tagged three-level patch-tree AMR with molecular transport. The MPI
work exponent is two so ownership reflects parabolic `r^2` subcycling cost.

Run it with, for example:

```sh
mpiexec -n 4 ./build-mpi/pelef_mpi_amr_reactive_1d \
  cases/mpi_sparse_amr_hotspot/hotspot.nml sparse_amr.csv
```
