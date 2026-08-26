"""
    scoring.jl

How a run is weighed against the Atlas LES: the per-variable normalization, the diagonal
observational noise, and the normalized misfit.

These are choices about an inversion rather than properties of the case — which variable
counts how much, and what uncertainty each carries — so they live here with the prior and
the observation map. `SOCRATES` owns reading the reference and aligning it with a run;
everything that assigns it a weight is below.
"""

function mean_nonzero_elements(x; all_zero = zero(eltype(x)))
    n = count(!iszero, x)
    n == 0 && return oftype(one(eltype(x)) * all_zero, all_zero)
    return sum(x) / n
end

"""
    ScoreTransform(; characteristic, obs_var_scaling, additional_uncertainty, uncertainty_floor)

Per-variable weights that put every compared quantity on a comparable scale.

`characteristic` is the magnitude a variable falls back to where the reference is all
zeros, `obs_var_scaling` the relative confidence in each variable, and the last two set the
diagonal noise the inversion adds.
"""
struct ScoreTransform{
    C <: AbstractDict, S <: AbstractDict, U <: AbstractDict, F <: AbstractDict,
}
    characteristic::C
    obs_var_scaling::S
    additional_uncertainty::U
    uncertainty_floor::F
end

const DEFAULT_CHARACTERISTIC = Dict{String, Float64}(
    "clw" => 1.0e-4,
    "cli" => 1.0e-8,
    "husra" => 1.0e-7,
    "hussn" => 1.0e-7,
    "lwp" => 1.0e-1,
    "iwp" => 1.0e-5,
    "rwp" => 3.0e-4,
    "swp" => 2.0e-4,
)

const DEFAULT_OBS_VAR_SCALING = Dict{String, Float64}(
    "clw" => 1.0,
    "cli" => (1 / 2.5)^2,
    "husra" => 1.0,
    "hussn" => (1 / 2.0)^2,
    "lwp" => (1 / 0.5)^2,
    "iwp" => (1 / 3.0)^2,
    "rwp" => 1.0,
    "swp" => 1.0,
)

const DEFAULT_ADDITIONAL_UNCERTAINTY =
    Dict{String, Float64}(name => 0.1 for name in keys(DEFAULT_CHARACTERISTIC))

const DEFAULT_UNCERTAINTY_FLOOR =
    Dict{String, Float64}(name => 0.05 for name in keys(DEFAULT_CHARACTERISTIC))

ScoreTransform(;
    characteristic = DEFAULT_CHARACTERISTIC,
    obs_var_scaling = DEFAULT_OBS_VAR_SCALING,
    additional_uncertainty = DEFAULT_ADDITIONAL_UNCERTAINTY,
    uncertainty_floor = DEFAULT_UNCERTAINTY_FLOOR,
) = ScoreTransform(
    characteristic,
    obs_var_scaling,
    additional_uncertainty,
    uncertainty_floor,
)

_entry(d, name, what) =
    get(d, name) do
        error("ScoreTransform has no $what for variable `$name`.")
    end

"""
    pool_var(transform, name, reference)

The variance `name` is normalized by: the scaled mean of the nonzero reference values,
floored at the variable's characteristic magnitude.
"""
function pool_var(transform::ScoreTransform, name::AbstractString, reference)
    scaling = _entry(transform.obs_var_scaling, name, "obs_var_scaling")
    FT = eltype(reference)
    cv = FT(_entry(transform.characteristic, name, "characteristic value"))
    μ = mean_nonzero_elements(reference; all_zero = cv)
    value = max(scaling * μ^2, scaling * cv^2)
    return iszero(value) ? eps(FT) : FT(value)
end

normalized_characteristic(transform::ScoreTransform, name::AbstractString) =
    1 / sqrt(_entry(transform.obs_var_scaling, name, "obs_var_scaling"))

"""
    uncertainty_diagonal(transform, name, series)

The diagonal observational-noise variance for each row of `series`.
"""
function uncertainty_diagonal(
    transform::ScoreTransform,
    name::AbstractString,
    series::AbstractMatrix,
)
    factor = _entry(transform.additional_uncertainty, name, "additional_uncertainty")
    fraction = _entry(transform.uncertainty_floor, name, "uncertainty_floor")
    floor = fraction * normalized_characteristic(transform, name)
    magnitude =
        [SOCRATES.nanmean(abs.(view(series, i, :))) for i in axes(series, 1)]
    return @. (factor * magnitude)^2 + floor^2
end

"""
    compare_to_les(case, output_dir; transform, vars, window, bounds, period, reduction)

The normalized misfit between a run's diagnostics and the Atlas LES reference, per
compared variable.
"""
function compare_to_les(
    c::SOCRATES.SOCRATESCase,
    output_dir::AbstractString;
    transform::ScoreTransform = ScoreTransform(),
    vars = SOCRATES.SCORED_VARS,
    window = score_window(c),
    bounds = default_z_bounds(c),
    period::Union{String, Nothing} = nothing,
    reduction::Union{String, Nothing} = nothing,
)
    model = run_outputvars(output_dir, vars; period, reduction)
    reference = les_outputvars(c; vars)
    out = Dict{String, Any}()
    for name in vars
        levels = model_levels(model[name], bounds)
        m = windowed_time_mean(
            restrict_to_levels(model[name], levels), window,
        )
        r = windowed_time_mean(
            reference_on_levels(reference[name], levels), window,
        )
        mvals, rvals = vec(m.data), vec(r.data)
        length(mvals) == length(rvals) || error(
            "`$name` has $(length(mvals)) model values and $(length(rvals)) \
             reference values after resampling.",
        )
        pv = pool_var(transform, name, rvals)
        mn, rn = mvals ./ sqrt(pv), rvals ./ sqrt(pv)
        d = filter(isfinite, mn .- rn)
        mse = isempty(d) ? NaN : sum(abs2, d) / length(d)
        out[name] = (; mse, rmse = sqrt(mse), model = mn, reference = rn, pool_var = pv)
    end
    return out
end

function print_comparison(comparison; io = stdout)
    println(
        io, rpad("variable", 10), rpad("n", 5), rpad("rmse (norm)", 14),
        rpad("mean model", 14), rpad("mean ref", 14), "pool_var",
    )
    for name in SOCRATES.SCORED_VARS
        haskey(comparison, name) || continue
        c = comparison[name]
        mm = Statistics.mean(filter(isfinite, c.model))
        mr = Statistics.mean(filter(isfinite, c.reference))
        println(
            io, rpad(name, 10), rpad(length(c.model), 5),
            rpad(round(c.rmse; sigdigits = 4), 14),
            rpad(round(mm; sigdigits = 4), 14),
            rpad(round(mr; sigdigits = 4), 14),
            round(c.pool_var; sigdigits = 4),
        )
    end
    return nothing
end
