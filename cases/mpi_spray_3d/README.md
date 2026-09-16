# Periodic MPI mixed spray

Build with `PELEF_ENABLE_MPI=ON` and `PELEF_ENABLE_SPRAY=ON`, then run from a
fresh output directory (input paths may be absolute):

```sh
mpiexec -n 4 /path/to/build/pelef_mpi_spray_3d /path/to/cases/mpi_spray_3d/periodic_mixture.nml
```

This fixed H2/O2-mechanism application advances periodic x-slab gas and two-component
particles with P/2-C/2-T/2-H-T/2-C/2-P/2 splitting. The example uses synthetic
water/peroxide liquid properties; it is an implementation exercise, not a qualified
physical spray model. Global cells must provide at least two x cells per rank and
three y/z cells. Rank zero reads and broadcasts the `mpi_spray` namelist.

The default initial gas is 800 K, density 1 kg/m3 and mass fractions H2=0.03,
O2=0.23, N2=0.74. `gas_fractions(10)` follows the pinned mechanism order
H2,H,O,O2,OH,H2O,HO2,H2O2,AR,N2. Liquid species, densities, heat capacities,
boiling temperatures (at 101325 Pa), latent heats and injection mass fractions
are configurable pairs. Initial cone injection is performed once; there is no
scheduled injection yet. `number=0,total_mass=0` disables injection.

Time controls are `final_time`, `dt_max`, `cfl`, `max_steps`, `rtol`, `atol` and
`chemistry` (`implicit`, `explicit`, or `cvode` in a SUNDIALS build). The code
reduces hydro/transport time-step limits across ranks and retries rejected complete
updates with halved steps. `transport_flags(4)` selects viscosity, thermal
conduction, species diffusion and barodiffusion. Solver/limiter defaults are HLLC/MC.

Outputs use `output_prefix`:

- `.history.csv`: step, time, global parcel count, particle substeps, total
  gas-plus-liquid mass/momentum/energy and relative mass/energy drift. Comment
  headers retain the mechanism SHA and global ID high-water mark.
- `.rankNNNNNN.gas.csv`: local cells with global integer indices, temperature and
  conserved gas fields.
- `.rankNNNNNN.particles.csv`: IDs, position, velocity, per-drop mass/temperature,
  multiplicity and per-drop component masses in configured liquid order.

Existing files are not overwritten. CSV output is not a restart checkpoint.
Selected mechanisms, scheduled/multiport injection, distributed breakup/wall films
and SGS remain to be integrated in this executable.

## Save and resume

Set `checkpoint_path='step.chk'` to save the final accepted state. Set
`restart_path='step.chk'` in a subsequent run, choose a fresh `output_prefix`, and
increase `final_time` to continue. Initial injection is skipped on restart. Gas,
particles, component masses, budgets, time, step count and the ID high-water mark
are restored. The MPI rank count may change if the new partition still has at
least two x cells per rank. The implementation gathers checkpoint data on root;
it is intended for bounded cases, not production parallel I/O.

The checkpoint binds the mechanism, mixture evolution, geometry, liquid properties,
injection definition and numerical controls. These must match. Rank count,
checkpoint/output paths, `final_time` and `max_steps` may change. Invalid or
incompatible checkpoints are rejected before output history is created. Existing
checkpoint destinations are refused. Focused same-rank continuation matches final
CSV values exactly; changed-rank continuation is compared with floating-point
tolerance, not claimed bitwise identical.
