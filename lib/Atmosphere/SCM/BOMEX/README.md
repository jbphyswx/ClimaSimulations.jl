# BOMEX

Trade-wind cumulus over the tropical Atlantic, as a single column.

```julia
using BOMEX
BOMEX.run_case(BOMEX.case(); output_dir = "runs/bomex")
```

## Where the case lives

`cases/bomex.toml` carries every number, each with its source. `src/` builds the grid, the
`AtmosModel` and the `AtmosSimulation` from it, and reads a run's own output.

The vertical profiles are `AtmosphericProfilesLibrary`'s, named rather than restated: `θ_li`,
`q_tot`, `u` and prescribed `tke` for the initial state, and subsidence, `dTdt`, `dqtdt` and the
geostrophic wind for the forcing. Pressure is integrated hydrostatically from the case's surface
pressure against `θ_li` and `q_tot`.

`BOMEXCase` implements the `ClimaAtmos.Setups` generics — `center_initial_condition`,
`surface_condition`, `subsidence_forcing`, `large_scale_advection_forcing` and
`coriolis_forcing`.

## Parameters

`bomex_params` Layer your own on top:

```julia
BOMEX.bomex_params(Float64; params = "my_overrides.toml")
BOMEX.bomex_params(Float64; params = Dict("some_parameter" => Dict("value" => 1.0)))
```

The default toml is at `configs/BOMEX.toml`

## Forcing and radiation

The surface fluxes are prescribed, not computed from a bulk formula: `MoninObukhov` is given the
case's `θ` flux, `q` flux and friction velocity directly.

There is no radiation. BOMEX prescribes its radiative cooling as part of the large-scale
temperature tendency, so an interactive scheme would count it twice — which is why no radiative
flux appears in `default_diagnostic_vars`.

## Tests

`test/runtests.jl` is offline: the grid, the profiles over the whole column, and that each
`Setups` generic carries the case's own values rather than ClimaAtmos's `nothing` default.

`test/integration.jl` constructs and steps a column and is a separately-invoked job, absent from
`runtests.jl` — construction alone takes minutes and gigabytes, so it is not suite work.
