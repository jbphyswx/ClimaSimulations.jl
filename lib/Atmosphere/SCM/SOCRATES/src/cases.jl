"""
    cases.jl

SOCRATES case definitions, vertical grids, flight constants, and parameter composition.
"""

using Dates: Dates
using NCDatasets: NCDatasets as NC
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF
using ClimaAtmos: ClimaAtmos as CA
using ClimaComms: ClimaComms
using ClimaParams: ClimaParams as CP
using TOML: TOML

# --- SocratesCase Definition & Registry ----------------------------------------------------- #

"""
    SocratesCase(flight_number, forcing_type)

One Atlas LES case: SOCRATES research flight `flight_number` forced by `forcing_type`
(`SSCF.ObsForcing()` or `SSCF.ERA5Forcing()`).
"""
struct SocratesCase{FT_TYPE <: SSCF.AbstractForcingType}
    flight_number::Int
    forcing_type::FT_TYPE
end

"""
    forcing_type(sym::Symbol)

The `SSCF.AbstractForcingType` for `:Obs` or `:ERA5`.
"""
forcing_type(::Val{:Obs}) = SSCF.ObsForcing()
forcing_type(::Val{:ERA5}) = SSCF.ERA5Forcing()
forcing_type(sym::Symbol) =
    sym in (:Obs, :ERA5) ? forcing_type(Val(sym)) :
    error("Unknown SOCRATES forcing type `:$sym`; expected `:Obs` or `:ERA5`")
forcing_type(ft::SSCF.AbstractForcingType) = ft
forcing_type(case::SocratesCase) = case.forcing_type

"""Short forcing label, `:Obs` or `:ERA5` (`SSCF.symbol`)."""
forcing_label(case::SocratesCase) = SSCF.symbol(case.forcing_type)

"""
    case_name(case)

The canonical case name, e.g. `"RF09_Obs"` or `"RF09_ERA5"`. Used for output subdirectories
and observation identifiers.
"""
case_name(case::SocratesCase) =
    string("RF", lpad(case.flight_number, 2, '0'), "_", forcing_label(case))

Base.show(io::IO, case::SocratesCase) = print(io, "SocratesCase(", case_name(case), ")")

"""
    case(name)

Parse a case name such as `"RF09_Obs"` or `"RF09_ERA5"` into a [`SocratesCase`](@ref).
"""
function case(name::AbstractString)
    m = match(r"^RF(\d{1,2})_(Obs|ERA5)$", name)
    isnothing(m) && error(
        "Cannot parse SOCRATES case name `$name`; expected e.g. \"RF09_Obs\" or \"RF09_ERA5\"",
    )
    return validate(SocratesCase(parse(Int, m[1]), forcing_type(Symbol(m[2]))))
end

case(c::SocratesCase) = c
case(flight_number::Integer, ft) =
    validate(SocratesCase(Int(flight_number), forcing_type(ft)))

# Backward-compatible alias
const socrates_case = case

"""
    all_cases()

Every valid (flight, forcing) Atlas LES case: the 5 Obs cases and the 6 ERA5 cases. Flight 11
has no Obs artifact, so it appears only under ERA5.
"""
function all_cases()
    cases = SocratesCase[]
    for ft in SSCF.forcing_types, flight in SSCF.flight_numbers
        SSCF.is_valid_flight_number(ft, flight) && push!(cases, SocratesCase(flight, ft))
    end
    return cases
end

# Backward-compatible alias
const socrates_cases = all_cases

"""
    validate(case)

Error unless SSCF has an Atlas artifact for this (flight, forcing) pair. Flight 11 has ERA5
forcing only.
"""
function validate(c::SocratesCase)
    SSCF.is_valid_flight_number(c.forcing_type, c.flight_number) || error(
        "No SOCRATES $(forcing_label(c)) artifact for flight $(c.flight_number) \
         (case $(case_name(c))). Valid Obs flights: \
         $(filter(f -> SSCF.is_valid_flight_number(SSCF.ObsForcing(), f), SSCF.flight_numbers)); \
         valid ERA5 flights: \
         $(filter(f -> SSCF.is_valid_flight_number(SSCF.ERA5Forcing(), f), SSCF.flight_numbers)).",
    )
    return c
end

# --- Run Length & Droplet Number Constants -------------------------------------------------- #

"""
Atlas LES run length [s] by forcing source: Obs cases run 12 h, ERA5 cases 14 h
(Atlas et al. 2020).
"""
const RUN_DURATION_SECONDS = Base.ImmutableDict(:Obs => 12 * 3600.0, :ERA5 => 14 * 3600.0)

"""
    t_end(case)

Simulation end time [s] for `case`, matching the Atlas LES run length.
"""
t_end(c::SocratesCase) = RUN_DURATION_SECONDS[forcing_label(c)]

"""
Prescribed cloud droplet number concentration [m^-3] by flight, from the Atlas LES `Nd` used
for each case. Enters the parameter set as `prescribed_cloud_droplet_number_concentration`.
"""
const N_CCN = Base.ImmutableDict(
    1 => 75.0e6,
    9 => 190.0e6,
    10 => 55.0e6,
    11 => 115.0e6,
    12 => 210.0e6,
    13 => 180.0e6,
)

"""
    n_ccn(case)

Prescribed cloud droplet number concentration [m^-3] for this case's flight.
"""
function n_ccn(c::SocratesCase)
    haskey(N_CCN, c.flight_number) ||
        error("No prescribed N_CCN for flight $(c.flight_number)")
    return N_CCN[c.flight_number]
end

# --- Clocks & Scoring Windows -------------------------------------------------------------- #

"""
    les_start_datetime(case)

True wall-clock start of the Atlas LES run (the Atlas reference time minus 12 h).
"""
les_start_datetime(c::SocratesCase) =
    SSCF.get_socrates_initial_time(c.flight_number)

"""
    simulation_start_date(case)

Epoch for the simulation clock. Fixed to 1970-01-01 so `t = 0` means "start of the run".
"""
simulation_start_date(::SocratesCase) = Dates.DateTime(1970, 1, 1)

"""Obs-case scoring window [s]: hours 10–12 of the 12 h run."""
const OBS_SCORE_WINDOW_SECONDS = (10 * 3600.0, 12 * 3600.0)

"""
    score_window(case)

The `(t_start, t_end)` window in seconds over which observations and model output are
time-averaged for this case.
"""
function score_window(c::SocratesCase)
    window =
        c.forcing_type isa SSCF.ObsForcing ? OBS_SCORE_WINDOW_SECONDS :
        era5_score_window(c.flight_number)
    return check_score_window(c, window)
end

"""
    check_score_window(case, (t0, t1))

Error unless the window is increasing and lies inside the case's run, `0 ≤ t0 < t1 ≤ t_end`.
"""
function check_score_window(c::SocratesCase, window)
    t0, t1 = window
    duration = t_end(c)
    (0 <= t0 < t1 <= duration) || error(
        "Scoring window ($t0, $t1) s for $(case_name(c)) is not within the run \
         [0, $duration] s (hours $(t0 / 3600) to $(t1 / 3600) of a \
         $(duration / 3600) h run).",
    )
    return (Float64(t0), Float64(t1))
end

"""
    era5_score_window(flight_number)

`(t_start, t_end)` in seconds for an ERA5 case, from the Atlas Table 2 metadata in
`SOCRATES_summary.nc`.
"""
function era5_score_window(flight_number::Integer)
    path = SSCF.atlas_socrates_summary_file(Int(flight_number))
    isfile(path) || error("Atlas SOCRATES summary file not found: $path")
    return NC.NCDataset(path, "r") do ds
        flights = vec(Array(ds["flight_number"]))
        i = findfirst(==(flight_number), flights)
        isnothing(i) && error(
            "Flight $flight_number not present in $path (has flights $(collect(flights)))",
        )
        bnds_var = ds["time_bnds"]
        size(bnds_var) == (2, length(flights)) || error(
            "$path `time_bnds` has size $(size(bnds_var)); expected \
             (2, $(length(flights))) as (nbnds, flight_number).",
        )
        bnds = bnds_var[:, i]
        reference = ds["reference_time"][i]
        offsets = (Float64(Dates.value(Dates.Second(b - reference))) for b in bnds)
        t0, t1 = (o + 12 * 3600.0 for o in offsets)
        (t0, t1)
    end
end

# --- Vertical Grids ------------------------------------------------------------------------ #

"""
    native_z(case)

The Atlas LES level centres [m] for `case`.
"""
native_z(c::SocratesCase) = collect(Float64, SSCF.default_new_z(c.flight_number))

"""
    native_faces(case)

The Atlas LES cell faces [m] for `case` — the model's default vertical grid.
"""
native_faces(c::SocratesCase) = faces_from_centers(native_z(c))

"""
    z_max_default(case)

Domain top [m] for `case`: the top face of the LES grid.
"""
z_max_default(c::SocratesCase) = last(native_faces(c))

"""
    faces_from_centers(centers; surface = 0.0)

Cell faces [m] for a grid whose centres are `centers`, from `f₁ = surface` and `fᵢ₊₁ = 2 cᵢ − fᵢ`.
"""
function faces_from_centers(centers::AbstractVector; surface::Real = 0.0)
    isempty(centers) && error("faces_from_centers needs at least one centre")
    faces = Vector{Float64}(undef, length(centers) + 1)
    faces[1] = Float64(surface)
    for (i, c) in enumerate(centers)
        faces[i + 1] = 2 * Float64(c) - faces[i]
    end
    issorted(faces) || error(
        "The given centres are not cell midpoints of a grid starting at $surface m: the implied faces \
         are not increasing. Pass `faces` directly for such a grid.",
    )
    return faces
end

"""
    centers_from_faces(faces)

Cell centres [m] of a grid with the given `faces`.
"""
centers_from_faces(faces::AbstractVector) =
    [(Float64(faces[i]) + Float64(faces[i + 1])) / 2 for i in 1:(length(faces) - 1)]

"""
    coarsen_faces_to_dz_min(faces, dz_min)

`faces` with interior faces dropped so every cell is at least `dz_min` thick, keeping the bottom and top.
"""
function coarsen_faces_to_dz_min(faces::AbstractVector, dz_min)
    fs = collect(Float64, faces)
    length(fs) >= 2 || error("A grid needs at least two faces, got $(length(fs))")
    isnothing(dz_min) && return fs
    minimum(diff(fs)) >= dz_min && return fs
    kept = [first(fs)]
    for f in fs[2:(end - 1)]
        (f - last(kept)) >= dz_min && push!(kept, f)
    end
    (last(fs) - last(kept)) >= dz_min ? push!(kept, last(fs)) : (kept[end] = last(fs))
    length(kept) >= 2 || error(
        "dz_min = $dz_min m leaves no cells in a column of depth $(last(fs) - first(fs)) m.",
    )
    return kept
end

"""
    socrates_grid(FT, case; faces, dz_min, context)

The column grid for `case`, built from explicit cell faces.
"""
function socrates_grid(
    ::Type{FT},
    c::SocratesCase;
    faces::AbstractVector = native_faces(c),
    dz_min = nothing,
    context = ClimaComms.context(),
) where {FT <: AbstractFloat}
    zf = coarsen_faces_to_dz_min(faces, dz_min)
    issorted(zf) ||
        error("Cell faces must be increasing; got $(length(zf)) faces spanning $(extrema(zf)) m.")
    domain = CA.CC.Domains.IntervalDomain(
        CA.CC.Geometry.ZPoint(FT(first(zf))),
        CA.CC.Geometry.ZPoint(FT(last(zf)));
        boundary_names = (:bottom, :top),
    )
    z_mesh = CA.CC.Meshes.IntervalMesh(domain, CA.CC.Geometry.ZPoint.(FT.(zf)))
    return CA.ColumnGrid(
        FT;
        context,
        z_elem = length(zf) - 1,
        z_max = FT(last(zf)),
        z_mesh,
    )
end

"""
    socrates_z(grid)

Centre-level heights [m] of `grid`, ascending.
"""
function socrates_z(grid)
    center_space = CA.get_spaces(grid).center_space
    z = CA.CC.Fields.coordinate_field(center_space).z
    return vec(Array(parent(z)))
end

"""
    socrates_z(FT, case; kwargs...)

Centre-level heights [m] for `case`, building the grid on the way.
"""
socrates_z(::Type{FT}, c::SocratesCase; kwargs...) where {FT <: AbstractFloat} =
    socrates_z(socrates_grid(FT, c; kwargs...))

# --- Parameter Composition ----------------------------------------------------------------- #

const TERMINAL_VELOCITY_SCALING_PARAMS = (;
    :cloud_liquid_terminal_velocity_scaling_factor => :liquid,
    :cloud_ice_terminal_velocity_scaling_factor => :ice,
    :rain_terminal_velocity_scaling_factor => :rain,
    :snow_terminal_velocity_scaling_factor => :snow,
)

terminal_velocity_scaling_defaults() = Dict{String, Any}(
    String(name) => Dict{String, Any}("value" => 1.0, "type" => "float") for
    name in keys(TERMINAL_VELOCITY_SCALING_PARAMS)
)

"""Base parameter TOML for the SOCRATES column."""
default_base_toml() = joinpath(dirname(@__DIR__), "configs", "base_prognostic_edmfx_1M.toml")

const ParamSource = Union{AbstractString, AbstractDict}

function _override_dict(source::AbstractString)
    isfile(source) || error("Parameter TOML not found: $source")
    endswith(source, ".toml") || error("Parameter source is not a .toml file: $source")
    return TOML.parsefile(source)
end
_override_dict(source::AbstractDict) =
    Dict{String, Any}(string(k) => v for (k, v) in source)

_source_list(params::ParamSource) = ParamSource[params]
_source_list(params::AbstractVector) = ParamSource[p for p in params]
_source_list(::Nothing) = ParamSource[]

"""
    n_ccn_override(case)

The case's prescribed cloud droplet number as a parameter-override dictionary.
"""
n_ccn_override(c::SocratesCase) = Dict{String, Any}(
    "prescribed_cloud_droplet_number_concentration" =>
        Dict{String, Any}("value" => n_ccn(c), "type" => "float"),
)

"""
    merge_param_sources(sources)

Merge parameter override `sources` left to right into one override dictionary.
"""
function merge_param_sources(sources)
    merged = Dict{String, Any}()
    for source in sources
        d = _override_dict(source)
        for (name, entry) in d
            merged[name] = entry
        end
    end
    return merged
end

"""
    socrates_toml_dict(FT, case; params, base_toml)

The `ClimaParams.ParamDict{FT}` for `case`: base defaults, terminal velocity defaults, case `n_ccn`,
then each entry of `params` in order.
"""
function socrates_toml_dict(
    ::Type{FT},
    c::SocratesCase;
    params = nothing,
    base_toml::AbstractString = default_base_toml(),
) where {FT <: AbstractFloat}
    sources = ParamSource[
        base_toml,
        terminal_velocity_scaling_defaults(),
        n_ccn_override(c),
        _source_list(params)...,
    ]
    return CP.create_toml_dict(FT; override_file = merge_param_sources(sources))
end

"""
    socrates_params(FT, case; params, base_toml, microphysics_model)

`ClimaAtmosParameters` for `case` in float type `FT`.
"""
socrates_params(
    ::Type{FT},
    c::SocratesCase;
    params = nothing,
    base_toml::AbstractString = default_base_toml(),
    microphysics_model = CA.NonEquilibriumMicrophysics1M(),
) where {FT <: AbstractFloat} = socrates_params(
    socrates_toml_dict(FT, c; params, base_toml);
    microphysics_model,
)

"""
    socrates_params(toml_dict; microphysics_model)

`ClimaAtmosParameters` from an already-composed `toml_dict`.
"""
socrates_params(
    toml_dict::CP.ParamDict;
    microphysics_model = CA.NonEquilibriumMicrophysics1M(),
) = CA.ClimaAtmosParameters(toml_dict; microphysics_model)
