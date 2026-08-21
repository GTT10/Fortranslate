# Decision 0017: explicit physical faces and mirrored wall states

The 0.18.0 path stores x fluxes on `0:nx` and y fluxes on `0:ny`. The exact
0.17.0 periodic hydro implementation remains the compatibility path. Physical
domains use typed ghost sampling, pressure-only impermeable inviscid wall
fluxes, mirrored wall velocity/temperature, and zero solid-wall species flux.
Fixed inflow is the configured initial mixture; outflow is constant
extrapolation. Catalytic walls and NSCBC remain outside this milestone.
