"""
    scoring.jl

Atlas LES reference profiles, the scored vertical region, and the normalized
misfit between a run and the reference.
"""

# --- Atlas LES reference ---------------------------------------------------- #

"""
The Atlas LES counterpart of each ClimaAtmos diagnostic, with the conversion from the
units [`atlas_var_specs`](@ref) produces onto the units ClimaAtmos writes.

Units and the mixing-ratio conversion are the registry's business, so this is a name
pairing plus whatever residual scaling the two vocabularies disagree on. `atlas_units` and
`model_units` are both pinned: a change to either side breaks `test/scoring.jl` rather
than silently rescaling an observation.

`lwp` pairs with `CWP` rather than `LWP`: `CWP` carries a 0.16% bias against the compared
quantity where `LWP` carries about 1%, and correlates at +0.84 with the precipitating area.

`hur` has no entry: ClimaAtmos declares no units for it, so the pairing against `RELH`
cannot be checked.
"""
const MODEL_TO_ATLAS = Dict{String, @NamedTuple{
    atlas::Symbol, atlas_units::String, model_units::String, to_model::Function,
}}(
    "ta" => (atlas = :TABS, atlas_units = "K", model_units = "K", to_model = identity),
    "ts" => (atlas = :SST, atlas_units = "K", model_units = "K", to_model = identity),
    "pfull" =>
        (atlas = :PRES, atlas_units = "Pa", model_units = "Pa", to_model = identity),
    "rhoa" => (
        atlas = :RHO, atlas_units = "kg/m3", model_units = "kg m^-3",
        to_model = identity,
    ),
    "hus" => (
        atlas = :QV, atlas_units = "kg/kg", model_units = "kg kg^-1",
        to_model = identity,
    ),
    "clw" => (
        atlas = :QCL, atlas_units = "kg/kg", model_units = "kg kg^-1",
        to_model = identity,
    ),
    "cli" => (
        atlas = :QCI, atlas_units = "kg/kg", model_units = "kg kg^-1",
        to_model = identity,
    ),
    "husra" => (
        atlas = :QPL, atlas_units = "kg/kg", model_units = "kg kg^-1",
        to_model = identity,
    ),
    "hussn" => (
        atlas = :QPI, atlas_units = "kg/kg", model_units = "kg kg^-1",
        to_model = identity,
    ),
    "lwp" => (
        atlas = :CWP, atlas_units = "kg/m2", model_units = "kg m^-2",
        to_model = identity,
    ),
    "iwp" => (
        atlas = :IWP, atlas_units = "kg/m2", model_units = "kg m^-2",
        to_model = identity,
    ),
    "rwp" => (
        atlas = :RWP, atlas_units = "kg/m2", model_units = "kg m^-2",
        to_model = identity,
    ),
    "swp" => (
        atlas = :SWP, atlas_units = "kg/m2", model_units = "kg m^-2",
        to_model = identity,
    ),
    "ua" =>
        (atlas = :U, atlas_units = "m/s", model_units = "m s^-1", to_model = identity),
    "va" =>
        (atlas = :V, atlas_units = "m/s", model_units = "m s^-1", to_model = identity),
    "tke" => (
        atlas = :TKE, atlas_units = "m2/s2", model_units = "m^2 s^-2",
        to_model = identity,
    ),
    "cdnc" => (
        atlas = :NCMN, atlas_units = "1/m3", model_units = "m^-3", to_model = identity,
    ),
    # ClimaAtmos writes cloud fraction as a percentage where the archive is a fraction
    "cl" => (atlas = :CLD, atlas_units = "1", model_units = "%", to_model = x -> 100 .* x),
    "hfls" => (
        atlas = :LHF, atlas_units = "W/m2", model_units = "W m^-2", to_model = identity,
    ),
    "hfss" => (
        atlas = :SHF, atlas_units = "W/m2", model_units = "W m^-2", to_model = identity,
    ),
    "pr" => (
        atlas = :PREC, atlas_units = "kg/m2/s", model_units = "kg m^-2 s^-1",
        to_model = identity,
    ),
    "rld" => (
        atlas = :RADLWDN, atlas_units = "W/m2", model_units = "W m^-2",
        to_model = identity,
    ),
    "rlu" => (
        atlas = :RADLWUP, atlas_units = "W/m2", model_units = "W m^-2",
        to_model = identity,
    ),
    "rsd" => (
        atlas = :RADSWDN, atlas_units = "W/m2", model_units = "W m^-2",
        to_model = identity,
    ),
    "rsu" => (
        atlas = :RADSWUP, atlas_units = "W/m2", model_units = "W m^-2",
        to_model = identity,
    ),
    # no downwelling longwave at the top, so the net there is the outgoing flux
    "rlut" => (
        atlas = :LWNTOA, atlas_units = "W/m2", model_units = "W m^-2",
        to_model = identity,
    ),
)

"""
    les_outputvars(case; vars, params)

The Atlas LES reference for `case`, as a `Dict` of ClimaAtmos short name to
`ClimaAnalysis.OutputVar`, in the units ClimaAtmos writes that name in.
"""
function les_outputvars(
    c::SOCRATES.SOCRATESCase;
    vars = SOCRATES.SCORED_VARS,
    params = SOCRATES.socrates_params(Float64, c),
)
    for name in vars
        haskey(MODEL_TO_ATLAS, name) || error(
            "`$name` has no Atlas counterpart in `MODEL_TO_ATLAS`, which pairs \
             $(join(sort(collect(keys(MODEL_TO_ATLAS))), ", ")).",
        )
    end
    raw = SOCRATES.read_atlas(c, [MODEL_TO_ATLAS[name].atlas for name in vars]; params)
    return Dict{String, ClimaAnalysis.OutputVar}(
        name => _outputvar(name, raw) for name in vars
    )
end

# A time series has one axis and a profile two, so the rank settles which dimensions the
# variable carries rather than a declaration that could disagree with the archive.
function _outputvar(name, raw)
    pair = MODEL_TO_ATLAS[name]
    values = pair.to_model(raw.data[pair.atlas])
    attribs = Dict{String, Any}(
        "short_name" => name,
        "units" => pair.model_units,
        "long_name" => "Atlas LES $name",
    )
    if ndims(values) == 2
        size(values) == (length(raw.z), length(raw.time)) || error(
            "Atlas variable `$(pair.atlas)` has size $(size(values)); expected \
             $((length(raw.z), length(raw.time))).",
        )
        dims = ClimaAnalysis.Var.OrderedDict("z" => raw.z, "time" => raw.time)
        dim_attribs = ClimaAnalysis.Var.OrderedDict(
            "z" => Dict("units" => "m"),
            "time" => Dict("units" => "s"),
        )
    else
        length(values) == length(raw.time) || error(
            "Atlas variable `$(pair.atlas)` has length $(length(values)); expected \
             $(length(raw.time)).",
        )
        dims = ClimaAnalysis.Var.OrderedDict("time" => raw.time)
        dim_attribs = ClimaAnalysis.Var.OrderedDict("time" => Dict("units" => "s"))
    end
    return ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, values)
end

# --- Resampling onto model levels ------------------------------------------- #

"""Relative tolerance on the edges of the reference column when resampling."""
const DEFAULT_EDGE_RTOL = 1.0e-4

"""
    reference_on_levels(var, levels; rtol)

`var` resampled in altitude onto `levels` [m].
"""
function reference_on_levels(
    var::ClimaAnalysis.OutputVar,
    levels;
    rtol::Real = DEFAULT_EDGE_RTOL,
)
    ClimaAnalysis.has_altitude(var) || return var
    padded = _pad_to_cell_extent(var)
    zp = padded.dims[ClimaAnalysis.altitude_name(padded)]
    lo, hi = first(zp), last(zp)
    tol_lo, tol_hi = rtol * (zp[2] - lo), rtol * (hi - zp[end - 1])
    requested = collect(Float64, levels)
    outside = filter(v -> v < lo - tol_lo || v > hi + tol_hi, requested)
    isempty(outside) ||
        error("Levels $outside m lie outside the Atlas LES column ($lo to $hi m).")
    resampled = ClimaAnalysis.resampled_as(padded; z = clamp.(requested, lo, hi))
    return _with_altitude(resampled, requested)
end

function _with_altitude(var::ClimaAnalysis.OutputVar, z)
    z_name = ClimaAnalysis.altitude_name(var)
    dims = ClimaAnalysis.Var.OrderedDict(
        name => name == z_name ? collect(Float64, z) : collect(d) for
        (name, d) in var.dims
    )
    return ClimaAnalysis.OutputVar(var.attributes, dims, var.dim_attributes, var.data)
end

# Extend the reference by half a cell at each end so that model levels sitting
# just outside the reference centres can still be resampled.
function _pad_to_cell_extent(var::ClimaAnalysis.OutputVar)
    z_name = ClimaAnalysis.altitude_name(var)
    z = collect(Float64, var.dims[z_name])
    faces = SOCRATES.faces_from_centers(z)
    axis = var.dim2index[z_name]
    n = size(var.data, axis)
    data = cat(
        selectdim(var.data, axis, 1:1),
        var.data,
        selectdim(var.data, axis, n:n);
        dims = axis,
    )
    dims = ClimaAnalysis.Var.OrderedDict(
        name => name == z_name ? vcat(first(faces), z, last(faces)) : collect(d)
        for (name, d) in var.dims
    )
    return ClimaAnalysis.OutputVar(var.attributes, dims, var.dim_attributes, data)
end

# --- Reading a run ---------------------------------------------------------- #

"""
    run_outputvars(output_dir, vars; period, reduction)

The scored diagnostics of a run, as `ClimaAnalysis.OutputVar`s.

`period` and `reduction` default to `nothing`, letting `ClimaAnalysis` discover
them; it errors if a run wrote more than one, which is the case where naming them
is actually required.
"""
function run_outputvars(
    output_dir::AbstractString,
    vars = SOCRATES.SCORED_VARS;
    period::Union{String, Nothing} = nothing,
    reduction::Union{String, Nothing} = nothing,
)
    simdir = ClimaAnalysis.SimDir(active_output_dir(output_dir))
    return Dict{String, ClimaAnalysis.OutputVar}(
        name => ClimaAnalysis.get(simdir; short_name = name, reduction, period)
        for name in vars
    )
end

active_output_dir(dir::AbstractString) =
    isdir(joinpath(dir, "output_active")) ? joinpath(dir, "output_active") : dir

model_levels(var::ClimaAnalysis.OutputVar, bounds) =
    ClimaAnalysis.has_altitude(var) ?
    scored_levels(var.dims[ClimaAnalysis.altitude_name(var)], bounds) : nothing

function restrict_to_levels(var::ClimaAnalysis.OutputVar, levels)
    (isnothing(levels) || !ClimaAnalysis.has_altitude(var)) && return var
    z_name = ClimaAnalysis.altitude_name(var)
    own = collect(Float64, var.dims[z_name])
    idx = indexin(collect(Float64, levels), own)
    any(isnothing, idx) &&
        error("Levels requested are not levels of variable $z_name.")
    axis = var.dim2index[z_name]
    dims = ClimaAnalysis.Var.OrderedDict(
        name => name == z_name ? own[idx] : collect(d) for (name, d) in var.dims
    )
    return ClimaAnalysis.OutputVar(
        var.attributes,
        dims,
        var.dim_attributes,
        copy(selectdim(var.data, axis, idx)),
    )
end

function windowed_time_mean(var::ClimaAnalysis.OutputVar, window)
    t = var.dims[ClimaAnalysis.time_name(var)]
    any(ti -> first(window) <= ti <= last(window), t) || error(
        "No times inside scoring window $window s (record spans $(extrema(t)) s).",
    )
    out = ClimaAnalysis.window(
        var,
        ClimaAnalysis.time_name(var);
        left = first(window),
        right = last(window),
    )
    mean_var = ClimaAnalysis.average_time(out)
    mean_var.data isa AbstractArray && return mean_var
    return ClimaAnalysis.OutputVar(
        mean_var.attributes,
        mean_var.dims,
        mean_var.dim_attributes,
        fill(mean_var.data),
    )
end

