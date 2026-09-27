# Reacting common-step MOL iteration — 2026-09-27

The MPI common integrator now permits chemistry in MOL mode. Its first
iteration forms the Euler predictor. Later iterations retain the previous
reacted candidate, evaluate its new spatial source, average that source with
the old source and reintegrate chemistry from the old state. This follows the
reacting loop in pinned `Source/Advance.cpp::do_mol_advance`, rather than
forming another Euler predictor each time. Nonreacting MOL still performs
one old/new correction. The particle predictor-builder restriction remains.

## Completed checks

GNU Fortran 13.3.0 / OpenMPI 4.1.6, existing Debug MPI/spray build,
tests/SUNDIALS disabled:

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_mpi_support --parallel 2 > build/mpi-mol-reaction-build.log 2>&1
env -u LD_LIBRARY_PATH -u LIBRARY_PATH /usr/bin/mpifort -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -I build/heptane-simple-mpi-debug/mpi_modules -I build/heptane-simple-mpi-debug/modules tests/unit/test_mpi_mol_reaction_3d.F90 build/heptane-simple-mpi-debug/libpelef_mpi_support.a build/heptane-simple-mpi-debug/libpelef_core.a -o build/mpi-mol-reaction-unit/check-short
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none /usr/bin/mpiexec -n 2 build/mpi-mol-reaction-unit/check-short > build/mpi-mol-reaction-unit/run-short.log 2>&1
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Final commands passed. One focused two-rank executable uses eight uniform
periodic cells owned by rank zero, 1200 K, 1 kg/m³, mass fractions H2=0.1,
O2=0.2, N2=0.699 and H=0.001, over `1e-10 s`. A homogeneous native reactor
calculation establishes nonzero chemical change. One and two MOL iterations
produce exactly 8 and 16 reactor calls. Species densities and temperature
match that homogeneous calculation with observed maximum error zero. The
chemical ledger closes. Deliberate failure on reactor call 9, after the first
iteration, rejects collectively and preserves gas, temperature and all ledgers.

The original radical-free 100 ns fixture failed the requirement for a
nontrivial reaction. A radical-seeded 100 ns diagnostic was then stopped
while still running to keep the focused test short; that interrupted run is
not counted as passed. The final 0.1 ns interval retains the same nontrivial
change and accuracy gates. Artifacts `run.log`, `run2.log` and `run-short.log`
preserve these distinctions. CTest registration is added; no full matrix,
selected application rebuild, C++ reacting comparison or spatially nonuniform
reacting MOL qualification was run. No Python source changed.

## Independent Jet failure

The retry-enabled Jet continuation ended with exit 1 and `Tree step retry
limit exhausted` at `2026-09-27T09:49:08.909616+00:00`, after
`1441.4868042636663 s`. Twelve half-step retries did not establish successful
continuation. Checkpoints through step 2400 remain intact. A diagnostic replay
from that checkpoint, using the exact same original binary, was launched in
`build/decane-jet-failure-trace`. Its GDB script captures accepted time/step
and next injection at the terminal error. The copied inputs, checkpoint hashes
and original binary identity are recorded in `provenance.json`; this replay
does not change acceptance tolerances or claim the failure has been repaired.
