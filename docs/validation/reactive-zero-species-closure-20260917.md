# Nonnegative chemistry closure for an absent final species — 2026-09-17

Parent `f7f1b3c` plus the chemistry correction. Validation used the accompanying
selected MPI spray frontend worktree; its sources/hashes and logs are retained
under `build/spray-selected-mpi-checks/` and documented in the selected MPI record.

The 2D chemistry kernel (also used plane-wise by 3D) writes the final species as
density minus the sum of the other species. If that species is absent, floating
summation can make it negative even with nonnegative integrated mass fractions.
The 21-species FFCM1_Red MPI water/methanol fixture exposed HE=-2.220446049250313e-16
at global cell (4,3,1), failing the unchanged strict nonnegative-species check.
The failing CSV remains under
`build/spray-selected-mpi-gnu-debug/selected-application-2/run-15lfj9ar/`.

When the final residual is negative by at most 32 epsilon times density, the
kernel now restores that species from its integrated abundance and places the
closure residual on the most abundant species. Larger negative residuals or a
negative corrected carrier reject the candidate. Positive-residual behavior is
unchanged; density and total energy are unchanged. This avoids clipping species
mass without compensating the closure. Test tolerances were not widened.

After the fix, the selected application checker passes at initial 2 and 4 ranks,
including schedule and same/changed-rank restart variants. The fixed-mechanism
application CTests pass 2/2 in 15.22 s, no failures/skips. GNU Fortran 13.3.0 Debug,
Open MPI 4.1.6. Provenance, contract, Python syntax and diff checks pass. No full
chemistry/mechanism matrix or physical validation is claimed.
