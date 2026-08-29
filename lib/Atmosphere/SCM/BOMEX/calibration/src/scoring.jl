"""
    scoring.jl

What is compared, over what window, on which levels.

BOMEX ships no observational reference — it is an LES intercomparison case — so the reference is
an argument, never a default. [`read_reference`](@ref) says what it may be.
"""

"""
Variables scored by default.

A starting point to edit, not the set that must be used; pass `vars` to override. No claim is
made here about which variables discriminate between parameter sets for this case — that has not
been measured.
"""
const default_calibration_vars = ("clw", "cl", "husra")

"""
    score_window(c)

The `(start, stop)` seconds compared: the last hour of the run.

Chosen here, not taken from a source. How much of a BOMEX run should be excluded as spin-up is a
question about the case, and it has not been settled by measurement in this package.
"""
score_window(c::BOMEX.BOMEXCase) = (BOMEX.t_end(c) - 3600.0, BOMEX.t_end(c))

"""The `(bottom, top)` metres compared: the whole column the case defines."""
default_z_bounds(::BOMEX.BOMEXCase) = (0.0, Float64(BOMEX.z_max()))

"""Levels of `z_grid` inside `bounds`."""
function scored_levels(z_grid, bounds)
    lo, hi = bounds
    levels = filter(z -> lo <= z <= hi, collect(Float64, z_grid))
    isempty(levels) && error("No model levels inside the compared region $bounds m.")
    return levels
end

"""
    read_reference(source) -> (; z, data)

The reference profiles to score against.

`source` may be a `Dict` mapping `"z"` and each variable name to a vector, or a path to a NetCDF
file carrying the same. Nothing is shipped: this case has no observational reference, so passing
one is the caller's job and there is no default to fall back on.
"""
function read_reference(d::AbstractDict)
    haskey(d, "z") || error("A reference needs a `\"z\"` entry giving its heights.")
    z = collect(Float64, d["z"])
    data = Dict{String, Vector{Float64}}(
        String(k) => collect(Float64, v) for (k, v) in d if String(k) != "z"
    )
    for (name, v) in data
        length(v) == length(z) || error(
            "Reference `$name` has $(length(v)) values against $(length(z)) heights.",
        )
    end
    return (; z, data)
end

function read_reference(path::AbstractString)
    isfile(path) || error("No reference file at $path.")
    return NC.NCDataset(path, "r") do ds
        haskey(ds, "z") || error("$path has no `z` variable.")
        z = collect(Float64, Array(ds["z"]))
        data = Dict{String, Vector{Float64}}()
        for name in keys(ds)
            name == "z" && continue
            v = Array(ds[name])
            ndims(v) == 1 && length(v) == length(z) &&
                (data[String(name)] = collect(Float64, v))
        end
        return (; z, data)
    end
end

"""
    reference_on_levels(z_src, values, levels)

`values` resampled onto `levels`, linear in height and refusing to extrapolate: a comparison
grid reaching past the reference's own column would otherwise be filled with invented values.
"""
function reference_on_levels(z_src, values, levels)
    z = collect(Float64, z_src)
    lo, hi = extrema(z)
    all(l -> lo - 1 <= l <= hi + 1, levels) || error(
        "Comparison levels $(extrema(levels)) m reach outside the reference column \
         $lo to $hi m.",
    )
    return CA.interp_vertical_prof(collect(Float64, levels), z, collect(Float64, values))
end

"""
    reference_values(source, names, levels)

Each of `names` from `source`, on `levels`.
"""
function reference_values(source, names, levels)
    ref = read_reference(source)
    out = Dict{String, Vector{Float64}}()
    for name in names
        haskey(ref.data, name) || error(
            "The reference has no `$name`; it carries " *
            join(sort(collect(keys(ref.data))), ", ") * ".",
        )
        out[String(name)] = reference_on_levels(ref.z, ref.data[name], levels)
    end
    return out
end

"""
    model_values(dir, names, levels; window, period)

Each of `names` from the run in `dir`, time-averaged over `window` and restricted to `levels`.

Both sides of a comparison are reduced the same way — mean over the window, then the same
levels — so a residual cannot come from averaging them differently.
"""
function model_values(
    dir::AbstractString,
    names,
    levels;
    window = nothing,
    period::Union{String, Nothing} = nothing,
)
    vars = BOMEX.run_outputvars(dir, collect(String, names); period)
    out = Dict{String, Vector{Float64}}()
    for name in names
        var = vars[String(name)]
        z = collect(Float64, var.dims[ClimaAnalysis.altitude_name(var)])
        t = collect(Float64, var.dims[ClimaAnalysis.time_name(var)])
        data = Array(var.data)
        # a run writes (time, z); put height first before reducing over time
        profile_by_time = size(data, 1) == length(t) ? permutedims(data, (2, 1)) : data
        keep =
            isnothing(window) ? eachindex(t) :
            findall(x -> first(window) <= x <= last(window), t)
        isempty(keep) && error(
            "No samples of `$name` inside $window s; the record spans $(extrema(t)) s.",
        )
        averaged = [Statistics.mean(view(profile_by_time, k, keep)) for k in axes(profile_by_time, 1)]
        idx = indexin(collect(Float64, levels), z)
        any(isnothing, idx) &&
            error("`$name` is not on the levels the observation was built on.")
        out[String(name)] = averaged[idx]
    end
    return out
end
