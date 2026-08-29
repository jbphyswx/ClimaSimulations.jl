"""
    BOMEXCase(FT; prognostic_tke, thermo_params)

The BOMEX case, carrying the vertical profiles it is initialised and forced from.

"""
struct BOMEXCase{FT <: AbstractFloat, P}
    prognostic_tke::Bool
    profiles::P
end

function BOMEXCase(
    ::Type{FT} = Float64;
    prognostic_tke::Bool = true,
    thermo_params = CA.TD.Parameters.ThermodynamicsParameters(FT),
) where {FT <: AbstractFloat}
    profiles = bomex_profiles(FT, thermo_params)
    return BOMEXCase{FT, typeof(profiles)}(prognostic_tke, profiles)
end

"""
    bomex_profiles(FT, thermo_params)

The case's height profiles: `θ_li`, `q_tot`, hydrostatic `p`, `u`, and prescribed `tke`.
"""
function bomex_profiles(::Type{FT}, thermo_params) where {FT <: AbstractFloat}
    θ = CA.APL.Bomex_θ_liq_ice(FT)
    q_tot = CA.APL.Bomex_q_tot(FT)
    p = CA.Setups.hydrostatic_pressure_profile(;
        thermo_params,
        p_0 = FT(reference_pressure()),
        θ,
        q_tot,
    )
    return (; θ, q_tot, p, u = CA.APL.Bomex_u(FT), tke = CA.APL.Bomex_tke_prescribed(FT))
end

"""The case's name, used to label a run's output directory."""
case_name(::BOMEXCase) = "BOMEX"

"""
    case(FT; kwargs...)

The BOMEX case. Takes the keyword arguments of [`BOMEXCase`](@ref).
"""
case(::Type{FT} = Float64; kwargs...) where {FT <: AbstractFloat} =
    BOMEXCase(FT; kwargs...)

"""Every case this package defines, which for BOMEX is the one."""
all_cases(::Type{FT} = Float64; kwargs...) where {FT <: AbstractFloat} =
    [case(FT; kwargs...)]

"""Run length [s] of `c`."""
t_end(::BOMEXCase) = t_end()

"""
    uniform_faces(top, n_cells)

The `n_cells + 1` cell faces [m] of a uniform column from the ground to `top`.
"""
function uniform_faces(top, n_cells)
    n_cells >= 1 || error("A column needs at least one cell, got $n_cells.")
    top > 0 || error("The domain top must be positive, got $top m.")
    return collect(range(0.0, Float64(top); length = Int(n_cells) + 1))
end

"""The cell faces [m] of the grid this package ships, from `cases/bomex.toml`."""
bomex_faces() = uniform_faces(z_max(), z_elem())

"""
    bomex_grid(FT, c; faces, context)

The column grid for `c`, built from explicit cell faces.
"""
function bomex_grid(
    ::Type{FT},
    ::BOMEXCase;
    faces::AbstractVector = bomex_faces(),
    context = ClimaComms.context(),
) where {FT <: AbstractFloat}
    zf = FT.(faces)
    issorted(zf) || error(
        "Cell faces must be increasing; got $(length(zf)) faces spanning $(extrema(zf)) m.",
    )
    domain = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint(FT(first(zf))),
        CC.Geometry.ZPoint(FT(last(zf)));
        boundary_names = (:bottom, :top),
    )
    z_mesh = CC.Meshes.IntervalMesh(domain, CC.Geometry.ZPoint.(zf))
    return CA.ColumnGrid(FT; context, z_elem = length(zf) - 1, z_max = FT(last(zf)), z_mesh)
end

"""Centre-level heights [m] of `grid`, ascending."""
function bomex_z(grid)
    center_space = CA.get_spaces(grid).center_space
    return vec(Array(parent(CC.Fields.coordinate_field(center_space).z)))
end

# --- Parameters ------------------------------------------------------------- #

"""
    reference_parameter_file()

The package's own ClimaParams overrides, `configs/BOMEX.toml`.

"""
reference_parameter_file() = joinpath(dirname(@__DIR__), "configs", "BOMEX.toml")

_override_dict(path::AbstractString) = TOML.parsefile(path)
_override_dict(d::AbstractDict) = Dict{String, Any}(d)
_override_dict(sources::AbstractVector) =
    reduce(merge!, map(_override_dict, sources); init = Dict{String, Any}())

"""
    bomex_toml_dict(FT; params)

The `ClimaParams.ParamDict{FT}` for the case: [`reference_parameter_file`](@ref), then each
entry of `params` in order, later sources winning.

`params` accepts a TOML path, an override dictionary, or a vector mixing both.
"""
function bomex_toml_dict(::Type{FT} = Float64; params = nothing) where {FT <: AbstractFloat}
    overrides = _override_dict(reference_parameter_file())
    isnothing(params) || merge!(overrides, _override_dict(params))
    return CP.create_toml_dict(FT; override_file = overrides)
end

"""`ClimaAtmosParameters` for the case."""
bomex_params(
    ::Type{FT} = Float64;
    params = nothing,
    microphysics_model = CA.NonEquilibriumMicrophysics1M(),
) where {FT <: AbstractFloat} =
    CA.ClimaAtmosParameters(bomex_toml_dict(FT; params); microphysics_model)