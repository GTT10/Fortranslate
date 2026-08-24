# Sparse MPI AMR reactive hotspot

This input drives the public `pelef_mpi_amr_reactive_1d` application through
solution-tagged three-level patch-tree AMR with molecular transport. The MPI
work exponent is two so ownership reflects parabolic `r^2` subcycling cost.

Run it with, for example:

```sh
mpiexec -n 4 ./build-mpi/pelef_mpi_amr_reactive_1d \
  cases/mpi_sparse_amr_hotspot/hotspot.nml sparse_amr.csv
```

## Checkpoint and rank-independent restart

`checkpoint_interval` selects a coarse-step cadence. A nonzero value requires
`checkpoint_file`. Setting `checkpoint_stop_after_write = .true.` makes the
driver exit cleanly after the first scheduled write. `restart_file` restores
the patch hierarchy, cell states, temperatures, physical time, step counters,
and AMR accounting before MPI ownership is rebuilt for the current rank count.

The checkpoint deliberately does not store the old rank-owner map, so a file
written by two ranks can be resumed by four or eight ranks:

```sh
mpiexec -n 2 ./build-mpi/pelef_mpi_amr_reactive_1d \
  cases/mpi_sparse_amr_hotspot/checkpoint_stop.nml
mpiexec -n 4 ./build-mpi/pelef_mpi_amr_reactive_1d \
  cases/mpi_sparse_amr_hotspot/restart.nml restarted.csv
```

The current format is a self-describing, portable formatted file gathered and
written by rank zero. It provides deterministic restart and rank
redistribution; scalable parallel checkpoint I/O remains future work.
