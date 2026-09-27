# Common particle MOL staging: qualification failure — 2026-09-27

The common particle wrapper now accepts `mol_mode`. It evaluates the old
kick/drift at stage zero and the new kick at stage one, retaining the new
particle source on later chemical iterations. This avoids both skipping the
new kick in nonreacting MOL and repeatedly moving particles during reaction
iterations. The MPI particle application forwards common MOL; recursive
particle MOL still rejects explicitly. Existing CTU/SDC staging is unchanged.

Pinned `Source/Advance.cpp::do_mol_advance` constructs old/new particle sources
once before later reaction iterations. `Source/Sources.cpp` routes the new
source to `particleMK`; `Source/Particle.cpp` permits its update when
`sub_ncycle==0`. Retaining particle exchange in later Fortran chemical
iterations preserves the coupled ledger; reacting particle equivalence to the
pinned reference has not been qualified.

## Completed checks and failed application

GNU Fortran 13.3.0 / OpenMPI 4.1.6; existing selected heptane Simple-transport
Debug MPI/spray build, tests/SUNDIALS disabled. These commands passed:

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2 > build/mpi-mol-particle-build.log 2>&1
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
python3 -m py_compile build/mpi-mol-particle/run_checks.py
git diff --check
```

`python3 build/mpi-mol-particle/run_checks.py` failed on its **first** launch:
the heated-wall four-parcel fixture rejected its first `1e-7 s` step. None of
its later planned restart checks ran. Preserve `build/mpi-mol-particle/driver.log`
and `build/mpi-mol-particle-final/continuous/fortran.log`. This is not a passed
particle MOL application qualification.

GDB replay scripts/logs under `build/mpi-mol-particle-debug` localize rejection
to the new-stage MOL hydro source, then `mol_reconstruct_face_3d`. At axis 1,
upper face, center fuel mass fraction `3.1336849956116898e-18`, density
`0.87897776141758854 kg/m³`, sound speed `400.62827443836892 m/s`, reconstruction
produces fuel partial density `-1.0376301946496e-27 kg/m³`. `candidate.log`
records the complete center/slopes before rejection. The pinned MOL formulas
also reconstruct partial densities from separately limited characteristic
components and have no explicit nonnegative partial-density guard at this
point; the Fortran implementation rejects the negative result.

No tolerance was loosened and no negative species state was accepted. A
positivity treatment and its effect on reference fidelity remain unresolved.
The initial GDB setup command also failed due to a wrong working directory;
the subsequent absolute-path scripts ran to their intended breakpoints.

Particle MOL, recursive particle MOL, reacting particle parity and full-duration
Spray-EB remain incomplete. The independent Jet failure is not diagnosed by
this distinct MOL fixture.
