"""
    diagnostics.jl

The diagnostics a comparison against this reference needs and ClimaAtmos has no name
for.
"""

_total_liquid(state, cache, time) =
    _total_liquid(state, cache, time, cache.atmos.microphysics_model)
_total_liquid(_, _, _, model) =
    error("`ql_all` needs a microphysics with prognostic rain; got $model")
_total_liquid(
    state, _, _,
    ::Union{CA.NonEquilibriumMicrophysics1M, CA.NonEquilibriumMicrophysics2M},
) = @. lazy((state.c.ρq_lcl + state.c.ρq_rai) / state.c.ρ)

_total_ice(state, cache, time) =
    _total_ice(state, cache, time, cache.atmos.microphysics_model)
_total_ice(_, _, _, model) =
    error("`qi_all` needs a microphysics with prognostic snow; got $model")
_total_ice(
    state, _, _,
    ::Union{CA.NonEquilibriumMicrophysics1M, CA.NonEquilibriumMicrophysics2M},
) = @. lazy((state.c.ρq_icl + state.c.ρq_sno) / state.c.ρ)

"""
    register_condensate_totals!()

Register `ql_all` and `qi_all`, the total liquid and total ice mixing ratios.

ClimaAtmos names the four species separately — `clw`/`cli` for cloud and
`husra`/`hussn` for precipitating — and has no profile for either total. Its
`clwvi`/`clivi` are column integrals, and cloud-only besides: they compute
`∫ρ(q_liq + q_ice)` and `∫ρ q_ice` whatever their comments say about including
precipitation.

Idempotent, because `add_diagnostic_variable!` errors on a name already registered
and this runs on every load of the package.
"""
function register_condensate_totals!(; registry = CA.Diagnostics.ALL_DIAGNOSTICS)
    haskey(registry, "ql_all") || CA.Diagnostics.add_diagnostic_variable!(
        short_name = "ql_all",
        units = "kg kg^-1",
        long_name = "Mass Fraction of Total Liquid Water",
        comments = "Cloud liquid plus rain, per mass of moist air.",
        compute = _total_liquid,
    )
    haskey(registry, "qi_all") || CA.Diagnostics.add_diagnostic_variable!(
        short_name = "qi_all",
        units = "kg kg^-1",
        long_name = "Mass Fraction of Total Ice",
        comments = "Cloud ice plus snow, per mass of moist air.",
        compute = _total_ice,
    )
    return nothing
end