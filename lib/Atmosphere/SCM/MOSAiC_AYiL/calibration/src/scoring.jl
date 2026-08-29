
#=

    What a calibration compares, and over what

    The variables, the window, the column and the normalization are calibration
    choices, not properties of the case — so they live here rather than in the
    package. The package supplies the pieces: `climaatmos_field` gives the
    reference's counterpart of a ClimaAtmos diagnostic, `column_water_path`
    integrates a path, and `best_simulation_top` bounds the usable column.

=#

"""
The default calibration vars
"""
const default_calibration_vars = ("ql_all", "qi_all", "lwp", "iwp")


const default_derive_calibration_vars = ()

"""
Which profile each scored water path integrates.

A path is a vertical integral, so the model's own `lwp` is integrated over the model
levels and the reference's `clwp_bar` over the reference's 286, and differencing the
two mixes a physical difference with an integration one. Both sides are rebuilt from
the `ρq` profiles on one grid instead (`docs/design.md` section 11).
"""
const water_paths = Dict(
    "lwp" => "ql_all", "iwp" => "qi_all",
    "clwp" => "clw", "rwp" => "husra", "clip" => "cli", "swp" => "hussn",
)

const default_calibration_water_paths =
    Tuple(v for v in default_calibration_vars if haskey(water_paths, v))

"""
Characteristic magnitude of each scored variable: the scale its misfit is normalized
by when the reference is empty.

Only used when a day's reference is identically zero, which happens often here — most
Arctic days have no liquid at all, or no ice. Being the wrong order of magnitude
against the simulation swamps the misfit, so these are not free.
"""
const DEFAULT_CHARACTERISTIC = Dict{String, Float64}(
    "ql_all" => 1.0e-3,
    "qi_all" => 1.0e-8,
    "clw" => 1.0e-4,
    "cli" => 1.0e-8,
    "husra" => 1.0e-2,
    "hussn" => 1.0e-7,
    "hus" => 1.0e-3,
    "lwp" => 1.0e-1,
    "iwp" => 1.0e-5,
)

"""
Per-variable weight on a variable's observation-error variance: the square of the
factor its error is divided by, so below 1 tightens that variable's block of the
noise covariance and weights it more.

EKP takes `(samples, covariances, names)` and does the whitening itself
(`EnsembleKalmanProcesses/src/Observations.jl:270-291`), building a low-rank
covariance from samples with `tsvd_cov_from_samples` (`:167-221`); rank truncation,
inflation, localization and the scheduler are set where the EKP process is
constructed, not here.
"""
const DEFAULT_OBS_VAR_SCALING = Dict{String, Float64}(
    "ql_all" => 1.0,
    "qi_all" => 1.0,
    "clw" => 1.0,
    "cli" => (1 / 2.5)^2,
    "husra" => 1.0,
    "hussn" => (1 / 2.0)^2,
    "hus" => 1.0,
    "lwp" => (1 / 0.5)^2,
    "iwp" => (1 / 3.0)^2,
)



"""
A map applied to a variable's physical values before it is normalized and compared.
The pooled statistic is taken of the transformed values, so a variable's entry in
[`DEFAULT_CHARACTERISTIC`](@ref) is stated in the space it is scored in.

Condensate spans orders of magnitude and is exactly zero over most of the column, so
comparing it linearly scores the one cloudy layer and nothing else.
"""
abstract type AbstractTransform end

struct NoTransform <: AbstractTransform end
(::NoTransform)(x) = x
Base.inv(::NoTransform) = NoTransform()

struct Log10OffsetTransform{FT} <: AbstractTransform
    offset::FT
end
# not `x::T`: the model's values may be Float32 where the offset is Float64
(transform::Log10OffsetTransform)(x) = log10(x + transform.offset)
Base.inv(transform::Log10OffsetTransform) = Base.Fix2(-, transform.offset) ∘ exp10

"""
Which variables are scored in a transformed space, and by what.

A transform defines the observable, so a transformed variable's `y` and its
observation covariance are both in the transformed space — nothing is transformed on
one side only. The offset is what keeps the log finite where the value is exactly
zero, so it, and not a separate characteristic magnitude, is what floors a
transformed variable.
"""
const default_variable_transformations = Dict(
    "cli" => Log10OffsetTransform(1.0e-8),
    "hussn" => Log10OffsetTransform(1.0e-8),
    "qi_all" => Log10OffsetTransform(1.0e-8),
    #
    "clip" => Log10OffsetTransform(1.0e-8),
    "swp" => Log10OffsetTransform(1.0e-6),
    "iwp" => Log10OffsetTransform(1.0e-6),

    # "clw" => Log10OffsetTransform(1.0e-10),
    # "husra" => Log10OffsetTransform(1.0e-10),
    # "ql_all" => Log10OffsetTransform(1.0e-10),
    # #
    # "clwp" => Log10OffsetTransform(1.0e-10),
    # "rwp" => Log10OffsetTransform(1.0e-10),
    # "lwp" => Log10OffsetTransform(1.0e-10),
)


"""
The fraction of a profile's own magnitude taken as its observational standard
deviation.
"""
const DEFAULT_UNCERTAINTY_FRACTION = 0.1

"""
    observation_variance(name, reference; characteristic, scaling, fraction)

The observation-error variance for `name`, in the space the variable is scored in:
`σ²` for `σ = max(characteristic floor, fraction × the profile's mean magnitude)`,
times the variable's [`DEFAULT_OBS_VAR_SCALING`](@ref) weight.

A band set on the observation, per profile, rather than the reference's own variance
over time. The reference is nudged, area-averaged and near-steady in the mean in stratocumulus,
so its time variance understates the uncertainty; and a model-error term would widen
the band everywhere instead of tied to the specific cloud variable magnitudes.

The mean of the nonzero values.

This is the diagonal EKP is handed as `covariances`; nothing rescales `y`.
"""
function observation_variance(
    name,
    reference;
    characteristic = DEFAULT_CHARACTERISTIC,
    scaling = DEFAULT_OBS_VAR_SCALING,
    fraction = DEFAULT_UNCERTAINTY_FRACTION,
    transform = NoTransform(),
    dof::Integer = length(reference),
    profile_dof::Integer = dof,
)
    key = String(name)
    cv = get(characteristic, key) do
        error("No characteristic magnitude for `$key`.")
    end
    factor = get(scaling, key) do
        error("No observation-variance scaling for `$key`.")
    end
    magnitude = mean_nonzero(reference; all_zero = 0.0)
    σ_physical = max(cv, fraction * magnitude)
    σ = σ_physical * transform_slope(transform, magnitude)
    # a profile contributes one term to the sum of squares per level and a water path
    # contributes one, so a scalar's variance is divided by a profile's length: before
    # the importance weight, a path then counts as much as a whole profile, and the
    # weights do not have to absorb a factor that changes with resolution
    value = factor * σ^2 * (dof / profile_dof)
    return iszero(value) ? eps(typeof(float(value))) : value
end

"""
    transform_slope(transform, x)

`|d transform / d x|` at `x`, which carries a band on the physical quantity into the
space the variable is scored in.

The band is always a statement about the physical value, so it is stated once and
propagated, rather than each transform needing its own notion of width. Where a
`log10` variable's magnitude is far below its offset the slope is large, widening the
band to say the observation carries no information there — which is right, since the
stored value is then the offset and not the quantity.
"""
transform_slope(::NoTransform, x) = one(float(x))
transform_slope(transform::Log10OffsetTransform, x) =
    1 / ((x + transform.offset) * log(10))


"""
    score_window(case)

The `(start, stop)` seconds a comparison averages over: the last hour of the run.
The reference's record is 300 s to `runtime`, so this sits inside it and excludes the
spin-up from the interpolated initial state.
"""
score_window(c::MOSAiC_AYiL.MOSAiCAYiLCase) =
    (MOSAiC_AYiL.t_end(c) - 3600.0, MOSAiC_AYiL.t_end(c))

"""
    z_bounds(case; tops)

The `(bottom, top)` metres a comparison covers: the ground to the day's calibration
top, above which the reference's ice is not reproducible (`docs/design.md` §12).
"""
function z_bounds(c::MOSAiC_AYiL.MOSAiCAYiLCase; tops = default_calibration_tops)
    top = get(tops, c.date) do
        error("AYiL day $(c.date) is not one of the $(length(tops)) calibrated days.")
    end
    return (0.0, top)
end

"""Levels of `z_grid` inside `bounds`."""
scored_levels(z_grid, bounds) =
    filter(z -> first(bounds) <= z <= last(bounds), collect(Float64, z_grid))

"""
    reference_on_levels(field, levels)

A reference `(; z, time, data)` resampled onto `levels`, one column per time.

Linear in height, and it errors rather than extrapolating: a comparison grid reaching
above the reference's own column would otherwise be filled with invented values.
"""
function reference_on_levels(field, levels)
    z = collect(Float64, field.z)
    lo, hi = extrema(z)
    all(l -> lo - 1 <= l <= hi + 1, levels) ||
        error(
            "Levels $(extrema(levels)) m reach outside the reference column \
             $(lo) to $(hi) m.",
        )
    out = Matrix{Float64}(undef, length(levels), size(field.data, 2))
    for j in axes(field.data, 2)
        itp = MOSAiC_AYiL.Intp.extrapolate(
            MOSAiC_AYiL.Intp.interpolate(
                (z,), Float64.(view(field.data, :, j)),
                MOSAiC_AYiL.Intp.Gridded(MOSAiC_AYiL.Intp.Linear()),
            ),
            MOSAiC_AYiL.Intp.Flat(),
        )
        out[:, j] = itp.(levels)
    end
    return out
end

"""
    windowed_mean(data, time, window)

The mean of a `(z, time)` array over the samples inside `window`.

Each reference variable carries its own time axis — the `modbulkmicrostat3` group is
one sample shorter than the rest — so the window is applied to the axis that came
with the data rather than to a shared one.
"""
function windowed_mean(data::AbstractMatrix, time::AbstractVector, window)
    keep = findall(t -> first(window) <= t <= last(window), time)
    isempty(keep) && error(
        "No samples inside the scoring window $window s; the record spans \
         $(extrema(time)) s.",
    )
    return [Statistics.mean(view(data, k, keep)) for k in axes(data, 1)]
end

"""The mean magnitude of the nonzero elements of `x`, or `all_zero` when there are none."""
function mean_nonzero(x; all_zero = zero(eltype(x)))
    values = filter(v -> isfinite(v) && !iszero(v), x)
    return isempty(values) ? all_zero : Statistics.mean(abs.(values))
end


"""
    model_field(output_dir, short_name; period, reduction)

One diagnostic of a run as `(; z, time, data)` with `data` a `(z, time)` array, so
both sides of a comparison have the same shape.
"""
function model_field(
    output_dir::AbstractString,
    short_name::AbstractString;
    period::Union{String, Nothing} = nothing,
    reduction::Union{String, Nothing} = nothing,
)
    var = only(
        values(
            MOSAiC_AYiL.run_outputvars(
                output_dir, (String(short_name),); period, reduction,
            ),
        ),
    )
    time = collect(Float64, var.dims[ClimaAnalysis.time_name(var)])
    ClimaAnalysis.has_altitude(var) ||
        return (; z = Float64[], time, data = reshape(collect(var.data), 1, :))
    z = collect(Float64, var.dims[ClimaAnalysis.altitude_name(var)])
    return (; z, time, data = permutedims(reshape(var.data, length(time), length(z)), (2, 1)))
end

function _misfit(
    name,
    mvals,
    rvals,
    characteristic,
    scaling = DEFAULT_OBS_VAR_SCALING,
    transforms = default_variable_transformations,
    profile_dof::Integer = length(rvals),
)
    length(mvals) == length(rvals) || error(
        "`$name` has $(length(mvals)) model values and $(length(rvals)) reference \
         values after resampling.",
    )
    t = get(transforms, String(name), NoTransform())
    mt, rt = t.(mvals), t.(rvals)
    # the band is stated on the physical reference and propagated through `t`, so the
    # variance is in the same space as the values it divides
    pv = observation_variance(
        name, rvals; characteristic, scaling, transform = t,
        dof = length(rvals), profile_dof,
    )
    mn, rn = mt ./ sqrt(pv), rt ./ sqrt(pv)
    d = filter(isfinite, mn .- rn)
    mse = isempty(d) ? NaN : sum(abs2, d) / length(d)
    return (;
        mse, rmse = sqrt(mse), model = mn, reference = rn, observation_variance = pv,
        model_physical = mvals, reference_physical = rvals,
    )
end

# The comparison grid's own cell edges, so a path is integrated over the cells the
# profiles live in rather than over the reference's.
function _faces_for(levels, faces)
    kept = Float64[]
    for (i, z) in enumerate(levels)
        j = findfirst(>=(z), faces)
        isnothing(j) && error("Level $z m is above the reference's top face.")
        push!(kept, Float64(faces[max(j - 1, 1)]))
        i == length(levels) && push!(kept, Float64(faces[j]))
    end
    return kept
end

"""Split `vars` into the profiles to fetch and the paths to integrate from them."""
function _split_vars(vars)
    names = String.(collect(vars))
    paths = filter(v -> haskey(water_paths, v), names)
    profiles = filter(v -> !haskey(water_paths, v), names)
    # a path is integrated from its profile, so that profile is fetched whether or not
    # it is scored in its own right
    return profiles, paths, unique(vcat(profiles, [water_paths[p] for p in paths]))
end

"""
    calibrated_values(case, levels, fetch; vars, window)

The physical value of each scored variable on `levels`, from whichever side `fetch`
reads.

`fetch(name) -> (; time, data)` supplies one variable already on `levels`; the
profiles are time-averaged over `window` and the paths integrated over the cells
`levels` sit in. One function for both sides, so a comparison cannot end up averaging
or integrating them differently.
"""
function calibrated_values(
    c::MOSAiC_AYiL.MOSAiCAYiLCase,
    levels,
    fetch;
    vars = default_calibration_vars,
    window = score_window(c),
)
    profiles, paths, needed = _split_vars(vars)
    values = Dict{String, Vector{Float64}}()
    for name in needed
        f = fetch(name)
        values[name] = windowed_mean(f.data, f.time, window)
    end
    out = Dict{String, Vector{Float64}}(name => values[name] for name in profiles)
    if !isempty(paths)
        # the face above the topmost level, so its cell is whole: truncating at the
        # level itself, or a metre above it, clips that cell and loses that much of
        # every path
        faces = MOSAiC_AYiL.truncate_faces_to_top(
            MOSAiC_AYiL.native_faces(c),
            MOSAiC_AYiL.MOSAiC_AYiL_face_above_center(last(levels)),
        )
        path_faces = _faces_for(levels, faces)
        ρ = let f = fetch("rhoa")
            windowed_mean(f.data, f.time, window)
        end
        for path in paths
            q = values[water_paths[path]]
            out[path] = MOSAiC_AYiL.column_water_path(
                reshape(q, :, 1), reshape(ρ, :, 1), path_faces,
            )
        end
    end
    return out
end

"""Read the reference, resampled onto `levels`."""
reference_fetch(c::MOSAiC_AYiL.MOSAiCAYiLCase, levels) =
    function (name)
        les = MOSAiC_AYiL.climaatmos_field(name, c.date)
        return (; les.time, data = reference_on_levels(les, levels))
    end

"""Read a run's diagnostics, restricted to `levels`."""
model_fetch(output_dir::AbstractString, levels; period = nothing, reduction = nothing) =
    function (name)
        m = model_field(output_dir, name; period, reduction)
        keep = indexin(collect(Float64, levels), m.z)
        any(isnothing, keep) &&
            error("`$name` is not on the levels the observation was built on.")
        return (; m.time, data = m.data[keep, :])
    end

"""
    comparison_levels(output_dir, bounds; period, reduction)

The levels a comparison runs on: a run's own, inside `bounds`.
"""
comparison_levels(
    output_dir::AbstractString,
    bounds;
    period = nothing,
    reduction = nothing,
) = scored_levels(
    model_field(output_dir, "rhoa"; period, reduction).z, bounds,
)

"""
    compare_to_les(case, output_dir; kwargs...)

The normalized misfit between a run's diagnostics and the DALES reference, per scored
variable and per water path.

The comparison grid is the run's own levels inside `bounds`; the reference is
resampled onto them, and both sides are averaged over `window` on their own time
axes. Water paths are integrated on that grid from both sides' `ρq`, never taken from
a stored path.
"""
function compare_to_les(
    c::MOSAiC_AYiL.MOSAiCAYiLCase,
    output_dir::AbstractString;
    vars = default_calibration_vars,
    window = score_window(c),
    bounds = z_bounds(c),
    period::Union{String, Nothing} = nothing,
    reduction::Union{String, Nothing} = nothing,
    levels = comparison_levels(output_dir, bounds; period, reduction),
    characteristic = DEFAULT_CHARACTERISTIC,
    scaling = DEFAULT_OBS_VAR_SCALING,
    transforms = default_variable_transformations,
)
    model = calibrated_values(
        c, levels, model_fetch(output_dir, levels; period, reduction); vars, window,
    )
    reference = calibrated_values(c, levels, reference_fetch(c, levels); vars, window)
    profile_dof = length(levels)
    return Dict{String, Any}(
        name => _misfit(
            name, model[name], reference[name], characteristic, scaling, transforms,
            profile_dof,
        ) for name in keys(model)
    )
end

"""The misfit table, in the order the variables were scored."""
function print_comparison(comparison; io = stdout, vars = default_calibration_vars)
    println(
        io, rpad("variable", 10), rpad("n", 5), rpad("rmse (norm)", 14),
        rpad("mean model", 14), rpad("mean ref", 14), "obs var",
    )
    for name in vars
        haskey(comparison, String(name)) || continue
        s = comparison[String(name)]
        println(
            io, rpad(name, 10), rpad(length(s.model), 5),
            rpad(round(s.rmse; sigdigits = 4), 14),
            rpad(round(Statistics.mean(filter(isfinite, s.model_physical)); sigdigits = 4), 14),
            rpad(round(Statistics.mean(filter(isfinite, s.reference_physical)); sigdigits = 4), 14),
            round(s.observation_variance; sigdigits = 4),
        )
    end
    return nothing
end