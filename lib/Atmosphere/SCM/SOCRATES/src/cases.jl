"""
    cases.jl

The SOCRATES case type, its vertical grid, and its parameter set.
"""

# --- Case ------------------------------------------------------------------- #

"""
    SocratesCase(flight_number, forcing_type)

One Atlas LES case: SOCRATES research flight `flight_number` forced by
`forcing_type` (`SSCF.ObsForcing()` or `SSCF.ERA5Forcing()`).
"""
struct SocratesCase{F <: SSCF.AbstractForcingType}
    flight_number::Int
    forcing_type::F
end

forcing_type(::Val{:Obs}) = SSCF.ObsForcing()
forcing_type(::Val{:ERA5}) = SSCF.ERA5Forcing()
forcing_type(ft::SSCF.AbstractForcingType) = ft
forcing_type(sym::Symbol) =
    sym in (:Obs, :ERA5) ? forcing_type(Val(sym)) :
    error("Unknown SOCRATES forcing type `:$sym`; expected `:Obs` or `:ERA5`")

"""Short forcing label, `:Obs` or `:ERA5`."""
forcing_label(c::SocratesCase) = SSCF.symbol(c.forcing_type)

"""
    case_name(case)

The canonical case name, e.g. `"RF09_Obs"`, used for output subdirectories and
observation identifiers.
"""
case_name(c::SocratesCase) =
    string("RF", lpad(c.flight_number, 2, '0'), "_", forcing_label(c))

Base.show(io::IO, c::SocratesCase) = print(io, "SocratesCase(", case_name(c), ")")

"""
    case(name)

Parse a case name such as `"RF09_Obs"` into a [`SocratesCase`](@ref).
"""
function case(name::AbstractString)
    m = match(r"^RF(\d{1,2})_(Obs|ERA5)$", name)
    isnothing(m) && error(
        "Cannot parse SOCRATES case name `$name`; expected e.g. \"RF09_Obs\".",
    )
    return validate(SocratesCase(parse(Int, m[1]), forcing_type(Symbol(m[2]))))
end

case(c::SocratesCase) = c
case(flight_number::Integer, ft) =
    validate(SocratesCase(Int(flight_number), forcing_type(ft)))

"""
    all_cases()

Every valid (flight, forcing) Atlas LES case. Flight 11 has no Obs artifact, so
it appears only under ERA5.
"""
all_cases() = [
    SocratesCase(flight, ft) for ft in SSCF.forcing_types,
    flight in SSCF.flight_numbers if SSCF.is_valid_flight_number(ft, flight)
]

"""
    validate(case)

Error unless SSCF has an Atlas artifact for this (flight, forcing) pair.
"""
function validate(c::SocratesCase)
    SSCF.is_valid_flight_number(c.forcing_type, c.flight_number) || error(
        "No SOCRATES $(forcing_label(c)) artifact for flight $(c.flight_number) \
         (case $(case_name(c))).",
    )
    return c
end

# --- Clocks ----------------------------------------------------------------- #

"""Simulation end time [s], matching the Atlas LES run length."""
t_end(c::SocratesCase) = run_duration(forcing_label(c))

"""The wall-clock start of the Atlas LES run, 12 h before the case reference time."""
les_start_datetime(c::SocratesCase) =
    SSCF.get_socrates_initial_time(c.flight_number)

"""Prescribed cloud droplet number concentration [m^-3] for this case's flight."""
n_ccn(c::SocratesCase) = droplet_number(c.flight_number)

"""
    score_window(case)

The `(t_start, t_end)` window [s] over which model and reference are time-averaged.
"""
score_window(c::SocratesCase) = check_score_window(
    c,
    c.forcing_type isa SSCF.ObsForcing ? obs_score_window() :
    era5_score_window(c.flight_number),
)

function check_score_window(c::SocratesCase, window)
    t0, t1 = window
    duration = t_end(c)
    (0 <= t0 < t1 <= duration) || error(
        "Scoring window ($t0, $t1) s for $(case_name(c)) is not within the run \
         [0, $duration] s.",
    )
    return (Float64(t0), Float64(t1))
end

"""
    era5_score_window(flight_number)

`(t_start, t_end)` [s] for an ERA5 case, from the Atlas Table 2 metadata in
`SOCRATES_summary.nc`.
"""
function era5_score_window(flight_number::Integer)
    path = SSCF.atlas_socrates_summary_file(Int(flight_number))
    isfile(path) || error("Atlas SOCRATES summary file not found: $path")
    return NC.NCDataset(path, "r") do ds
        flights = vec(Array(ds["flight_number"]))
        i = findfirst(==(flight_number), flights)
        isnothing(i) && error(
            "Flight $flight_number not present in $path (has $(collect(flights))).",
        )
        bnds = ds["time_bnds"]
        size(bnds) == (2, length(flights)) || error(
            "$path `time_bnds` has size $(size(bnds)); expected \
             (2, $(length(flights))).",
        )
        reference = ds["reference_time"][i]
        offsets = (Float64(Dates.value(Dates.Second(b - reference))) for b in bnds[:, i])
        Tuple(o + 12 * 3600.0 for o in offsets)
    end
end

# --- Vertical grid ---------------------------------------------------------- #

"""The Atlas LES level centres [m] for `case`."""
native_z(c::SocratesCase) = collect(Float64, SSCF.default_new_z(c.flight_number))

"""The Atlas LES cell faces [m] for `case` — the model's default vertical grid."""
native_faces(c::SocratesCase) = faces_from_centers(native_z(c))

"""Domain top [m] for `case`: the top face of the LES grid."""
z_max_default(c::SocratesCase) = last(native_faces(c))

"""
    faces_from_centers(centers; surface = 0.0)

Cell faces [m] of the grid whose centres are `centers`, from `f₁ = surface` and
`fᵢ₊₁ = 2 cᵢ − fᵢ`.
"""
function faces_from_centers(centers::AbstractVector; surface::Real = 0.0)
    isempty(centers) && error("faces_from_centers needs at least one centre")
    faces = Vector{Float64}(undef, length(centers) + 1)
    faces[1] = Float64(surface)
    for (i, c) in enumerate(centers)
        faces[i + 1] = 2 * Float64(c) - faces[i]
    end
    issorted(faces) || error(
        "The given centres are not cell midpoints of a grid starting at \
         $surface m: the implied faces are not increasing.",
    )
    return faces
end

"""Cell centres [m] of a grid with the given `faces`."""
centers_from_faces(faces::AbstractVector) =
    [(Float64(faces[i]) + Float64(faces[i + 1])) / 2 for i in 1:(length(faces) - 1)]

"""
    coarsen_faces_to_dz_min(faces, dz_min)

`faces` with interior faces dropped so every cell is at least `dz_min` thick,
keeping the bottom and top.
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
        "dz_min = $dz_min m leaves no cells in a column of depth \
         $(last(fs) - first(fs)) m.",
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
    issorted(zf) || error(
        "Cell faces must be increasing; got $(length(zf)) faces spanning \
         $(extrema(zf)) m.",
    )
    domain = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint(FT(first(zf))),
        CC.Geometry.ZPoint(FT(last(zf)));
        boundary_names = (:bottom, :top),
    )
    z_mesh = CC.Meshes.IntervalMesh(domain, CC.Geometry.ZPoint.(FT.(zf)))
    return CA.ColumnGrid(FT; context, z_elem = length(zf) - 1, z_max = FT(last(zf)), z_mesh)
end

"""Centre-level heights [m] of `grid`, ascending."""
function socrates_z(grid)
    center_space = CA.get_spaces(grid).center_space
    return vec(Array(parent(CC.Fields.coordinate_field(center_space).z)))
end

socrates_z(::Type{FT}, c::SocratesCase; kwargs...) where {FT <: AbstractFloat} =
    socrates_z(socrates_grid(FT, c; kwargs...))

# --- Parameters ------------------------------------------------------------- #

"""One parameter override source: a path to a TOML file, or a parsed override dictionary."""
const ParamSource = Union{AbstractString, AbstractDict}

_override_dict(source::AbstractDict) =
    Dict{String, Any}(string(k) => v for (k, v) in source)
function _override_dict(source::AbstractString)
    isfile(source) || error("Parameter TOML not found: $source")
    endswith(source, ".toml") ||
        error("Parameter source is not a .toml file: $source")
    return TOML.parsefile(source)
end

_source_list(::Nothing) = ParamSource[]
_source_list(source::ParamSource) = ParamSource[source]
_source_list(sources::AbstractVector) = ParamSource[s for s in sources]

"""
    socrates_toml_dict(FT, case; params)

The `ClimaParams.ParamDict{FT}` for `case`: [`SOCRATES_PARAMETERS`](@ref), then the
case's prescribed droplet number, then each entry of `params` in order, later sources
winning. No file is read unless `params` names one.

`params` accepts a TOML path, an override dictionary, or a vector mixing both, and can
override any value of either layer beneath it.
"""
function socrates_toml_dict(
    ::Type{FT},
    c::SocratesCase;
    params = nothing,
) where {FT <: AbstractFloat}
    sources = ParamSource[
        parameter_overrides(),
        Dict{String, Any}(
            "prescribed_cloud_droplet_number_concentration" =>
                Dict{String, Any}("value" => n_ccn(c), "type" => "float"),
        ),
        _source_list(params)...,
    ]
    overrides = Dict{String, Any}()
    for source in sources
        merge!(overrides, _override_dict(source))
    end
    return CP.create_toml_dict(FT; override_file = overrides)
end

"""`ClimaAtmosParameters` for `case`."""
socrates_params(
    ::Type{FT},
    c::SocratesCase;
    params = nothing,
    microphysics_model = CA.NonEquilibriumMicrophysics1M(),
) where {FT <: AbstractFloat} =
    socrates_params(socrates_toml_dict(FT, c; params); microphysics_model)

socrates_params(
    toml_dict::CP.ParamDict;
    microphysics_model = CA.NonEquilibriumMicrophysics1M(),
) = CA.ClimaAtmosParameters(toml_dict; microphysics_model)
