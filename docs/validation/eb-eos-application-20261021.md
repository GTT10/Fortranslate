# EB-tree caloric convergence input, 2026-10-21

The ordinary EB-tree namelist exposes `pinned_eos_convergence`, default false,
through the existing NASA7 policy configuration. MPI additionally agrees the
normalized thermodynamic policy across ranks before initialization. Existing
patch checkpoint policy records already bind the precision setting; tree
checkpoints retain that protection through their patch records.

The focused application check passes. It uses the existing hot H2/O2/N2
CVODE case with the selected full mechanism, reaction-off comparison,
regridding, owner rebalancing, serial/MPI field identity and two-to-three-rank
checkpoint continuation. Added refusal checks cover a changed precision
setting on restart and mixed precision across an MPMD two-rank launch.
Both refusal checks produce no trajectory output, and the restart refusal
preserves checkpoint bytes. Serial, MPI and changed-rank continuation fields
and synchronized snapshots are byte-identical. Hydrogen reacted fraction is
`2.2952255206e-4`; the serial maximum normalized budget residual is
`7.31e-16` and chemical H/O/N residuals are below `5.51e-16`.

Artifacts: `build/eb-eos-application-20261021/`. The first configure attempt
failed before compilation because the environment selected Intel MPI for
GNU Fortran. The `gnu-openmpi/` retry configured successfully but exposed
a missing selected-runtime build dependency for `mpi_reactive_eos_policy_mod`.
That source is now included in the selected MPI EB-tree runtime.
The `dependency-fixed/` attempt compiled but failed the existing runtime
dependency audit: CMake had selected Intel `libmpi.so` through inherited
library search paths despite the explicit OpenMPI wrapper. All failed
attempts remain intact. The successful `isolated-openmpi/` attempt removes
oneAPI compiler/library/pkg-config search variables only in its subprocess
environment, checks CMake MPI library paths before compiling, and uses a
clean Debug build with explicit `/usr/bin/gfortran`,
`/usr/bin/mpifort.openmpi`, `/usr/bin/mpiexec.openmpi` and SUNDIALS 7.2.
The h2o2_full pinned-source provenance/regeneration preflight passed.
`record.json` in each attempt records source hashes and exact commands;
`resume.json` identifies the active job and completion requirements.
Python syntax, `git diff --check` and the project contract pass; their exact
commands and exits are recorded in `source-checks.json`. An initial contract
CLI invocation omitted its required `--source` argument; its usage-error log
is retained separately from the successful check.

The successful wrapper exits zero after 223.43 s: configure 2.67 s, clean
build 174.46 s and application check 46.23 s. Exactly one CTest passes, with
zero failures or skips: `regression_selected_eb_tree_converged_cvode_restart`.
Its detailed report is
`isolated-openmpi/registration/tests/selected_eb_tree_converged_cvode_restart/reacting_restart_4kmxh0nv/checks.json`.
It records five successful runs and four expected refusals (output policy,
checkpoint alias, changed restart precision and mixed-rank precision).
The compiler is GNU Fortran 13.3.0; OpenMPI library resolution and the
post-link MPI runtime audit pass. Source hashes still match the completed run
when the result is audited. This checks the selected full-mechanism gas
frontends; the fixed-mechanism and spray/LES frontends were not independently
built or exercised in this application gate.

This application gate does not establish original-PeleC parity, asymptotic
convergence, spray/LES validity, a full build/install matrix or physical
validation. The separate completed mesh evidence is in the
[three-grid record](eb-three-grid-20261021.md).
