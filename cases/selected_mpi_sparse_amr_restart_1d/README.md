# Selected sparse MPI AMR 1D restart

This case qualifies a configure-time selected full-H2/O2 mechanism in the
sparse MPI AMR patch-tree application. The checkpoint is written with one MPI
rank and restored with two and four ranks; MPI ownership is intentionally not
part of the persisted schema.

The fixed and selected inputs use the same species composition, implicit
chemistry policy, mesh, AMR criteria, and final time. Their composite CSV files
must therefore agree exactly.
