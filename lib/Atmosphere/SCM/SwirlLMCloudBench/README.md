# SwirlLMCloudBenchSim

Single-column ClimaAtmos runs of the Swirl-LM CloudBench sites: 500 sites × 4 months
(Jan/Apr/Jul/Oct) × 5 experiments (`amip`, `amip_p4k`, `amip_4xco2`, `amip_p4k_2xco2`,
`amip_p4k_4xco2`), forced as in Shen et al. (2022) and documented by Chammas et al. (2026).

[`SwirlLMCloudBench`](https://github.com/jbphyswx/SwirlLMCloudBench.jl) owns the dataset and
its ClimaAtmos extension: the sounding, the forcing, the initial condition, the insolation,
the surface condition and the ClimaParams overrides all come from there. This package is the
column assembled around them.

## Running a case

```julia
using SwirlLMCloudBench: Simulation as S
using SwirlLMCloudBenchSim: SwirlLMCloudBenchSim as CB

case = S.CloudBenchInstance(0, 1, :amip)          # site 0, January, AMIP
CB.run_case(case; output_dir = "output/site0_jan_amip")
```

`cloudbench_simulation(FT, case; output_dir, kwargs...)` returns the `AtmosSimulation`
without solving it. Every default is a keyword: `grid`, `params`, `setup`, `t_end`, `dt`,
`diagnostics`, `ode_config`, `jacobian`, `callback_kwargs`, and `model_kwargs` for any field
of `AtmosModel`.

The case's own `parameters.json` supplies its SST, zenith angle and TOA insolation, and its
experiment fixes the CO₂, so nothing about the scenario has to be passed by hand.

## Grid

The reference LES is 6000 m deep on 480 levels, which is affordable for one case and not for
an ensemble:

```julia
CB.cloudbench_grid(Float64)                  # the LES resolution, 480 cells
CB.cloudbench_grid(Float64; dz_min = 200.0)  # coarsened, 30 cells
```

## Reading the published output

`data.zarr` is `(z, x, y, t)` at 480 × 124 × 124 × 73, about 2.2 GB per four-dimensional
field, and covers the last simulated day at 1200 s. `reference` reads named fields and
returns a column: the horizontal mean of a four-dimensional array, or the array as stored
if the store already reduced it.

```julia
store = CB.reference_store(case)
CB.reference_axes(store)                              # z [m] and time [s]
CB.reference(case, ("lwp", "olr", "q_c"); store)      # column-mean if 4-D
```

## Tests

```bash
julia --project=test test/runtests.jl     # no network
julia --project=test test/reference.jl    # reads data.zarr over HTTPS
julia --project=test test/integration.jl  # builds a column and steps it; network, minutes
```

Only `runtests.jl` is the default suite. The other two are gated jobs: they need the network,
and `integration.jl` also needs several GB and minutes of wall clock.

## Microphysics

CloudBench runs single-moment warm rain with instantaneous condensation. ClimaAtmos carries
no equilibrium 1-moment model, so `cloudbench_model` defaults to
`NonEquilibriumMicrophysics1M`, which relaxes condensation over a finite timescale. That
timescale is a ClimaParams entry, so it is a parameter change rather than a code change.
