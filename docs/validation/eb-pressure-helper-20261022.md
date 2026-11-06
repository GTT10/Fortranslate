# Transactional EB pressure representation, 2026-10-22

The pressure candidate now calls `reconstruct_eb_pressure_state_3d` in
`eb_reactive_transfer_3d_mod` at all three former probe-local reconstruction
sites (fine halo, coarse output, covered coarse spatial representation).
The common routine accepts a target pressure and returns the numerical energy
density correction. Density, momentum and species remain bit-identical;
total energy and held internal energy receive the same correction. A failed
EOS admission publishes none of state, held energy, temperature or correction.
Arithmetic bounds precede potentially overflowing intermediate operations.
The routine is explicitly for nonrotating auxiliary states, not a conservative
update of valid leaf cells.

## Remaining application connection

Ordinary serial and MPI halo fills currently reconstruct gas and held energy
then call `recover_reactive_eb_cell_3d`. Connecting this helper requires a
separate pressure field formed from fluid-volume-weighted fine EOS pressures
on covered donors, including temporally interpolated ancestors. Both the
coarse covered spatial representation and fine halos must consume that same
field. MPI must transfer the pressure stencil with the gas/held stencil and
collectively refuse failures. Applying only a fine-halo pressure repair would
repeat the already rejected asymmetric candidate.

The pressure policy also needs held-stage/replay identity and checkpoint
admission before application exposure, and explicit accounting if any accepted
valid-leaf energy is changed. This extraction does not enable pressure
reconstruction in ordinary applications. The velocity option already has its
own application integration evidence.

## Evidence

Artifacts: `build/eb-pressure-helper-20261022/`.
Validation results and exact source/binary identity are recorded in `record.json`.

The focused `unit_eb_reactive_transfer_3d` CTest passes (1/1, zero failures or
skips). It checks the pressure target, immutable material/momentum bits, a
nonzero total-minus-held gap, returned correction, and complete rollback for
negative/NaN pressure, EOS temperature-range refusal, NaN held energy or
initial temperature, and zero density. Existing transfer tests also execute.
The pressure diagnostic Python tests pass (21/21).

Build command:

```sh
cmake --build build/eb-eos-application-20261021/isolated-openmpi/registration \
  --target test_eb_reactive_transfer_3d test_eb_tree_reacting_reference -j 4
ctest --test-dir build/eb-eos-application-20261021/isolated-openmpi/registration \
  -R '^unit_eb_reactive_transfer_3d$' --output-on-failure
python3 -m unittest discover -s tests/python -p '*pressure*.py'
```

The existing GNU Fortran Debug/OpenMPI build is reused without changing the
compiler or optional backends. The linker retains executable-stack warnings
for existing Fortran internal-procedure objects. This is a focused numerical
extraction check, not a new full application/MPI, install or physical gate.

Both frozen binaries completed the four-step 50 ns cut-cell controls for
smooth isothermal and imported nonuniform-temperature states, with pressure
and velocity policies. All 32 corresponding CSV artifacts are byte-identical,
including initial interface/ghost diagnostics, evolution, cell histories and
synchronization outputs. Both runs pass the existing budget, localization and
synchronization checks. This verifies extraction parity, not improvement of
the imported nonisothermal pressure residual.

The exact control commands and input/binary SHA256 values are nested in
`record.json`; source hashes are in `source-identity.json`. Provenance, Python
syntax, project contract and `git diff --check` passed. No tolerance or golden
hash was changed.
