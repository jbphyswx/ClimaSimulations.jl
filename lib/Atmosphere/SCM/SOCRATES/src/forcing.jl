"""
    forcing.jl

The SOCRATES large-scale forcing: a forcing type owned by this package, with
`external_forcing_cache` / `external_forcing_tendency!` methods composing
ClimaAtmos's shared forcing kernels.
"""

"""
SSCF forcing variables consumed by [`SOCRATESForcing`](@ref), in ClimaAtmos's
canonical terms: `dTdt_hadv`/`dqtdt_hadv` are the large-scale horizontal
advective tendencies (`tntha`/`tnhusha`), `T_nudge`/`qt_nudge`/`u_nudge`/`v_nudge`
the relaxation targets (`ta`/`hus`/`ua`/`va`), and `subsidence` the large-scale
vertical velocity (`wa`), positive upward.

SSCF carries no vertical eddy-fluctuation tendency, so there is no
`tntva`/`tnhusva` counterpart: vertical transport comes from `subsidence` acting
on the model's evolving profiles.
"""
const SSCF_FORCING_VARS =
    (:dTdt_hadv, :dqtdt_hadv, :T_nudge, :qt_nudge, :u_nudge, :v_nudge, :subsidence)

"""Surface series SSCF supplies, carried alongside the column profiles."""
const SSCF_SURFACE_VARS = (:pg, :Tg, :Tsfc, :qg, :qsfc)

"""
    SOCRATESForcing(FT, source; scalar_nudge_timescale, wind_nudge_timescale)

Large-scale forcing for one SOCRATES case: horizontal advection of temperature
and total water, relaxation of temperature, total water and horizontal wind
toward the case profiles, and large-scale subsidence.

`source` is either a [`SOCRATESCase`](@ref), sampled from SSCF, or the path of a
file written by [`write_forcing_file`](@ref) holding those same profiles. All
grid-dependent work happens once in `ClimaAtmos.external_forcing_cache`.

The relaxation rate is uniform in height. ClimaAtmos's `DefaultTimescale()` would
instead apply the GCM-driven ramp, which is zero below
`gcmdriven_relaxation_minimum_height` and would leave the boundary layer — the
thing being modelled — unconstrained.

# Keyword Arguments

  - `scalar_nudge_timescale`: Relaxation timescale of `T` and `q_tot` [s].
  - `wind_nudge_timescale`: Relaxation timescale of `u` and `v` [s].
"""
struct SOCRATESForcing{FT, S <: Union{SOCRATESCase, AbstractString}}
    source::S
    inv_τ_scalar::FT
    inv_τ_wind::FT
end

_validated_source(c::SOCRATESCase) = validate(c)
function _validated_source(path::AbstractString)
    isfile(path) || error("SOCRATES forcing file not found: $path")
    return path
end

function SOCRATESForcing(
    ::Type{FT},
    source;
    scalar_nudge_timescale,
    wind_nudge_timescale,
) where {FT <: AbstractFloat}
    checked = _validated_source(source)
    for (name, τ) in (
        ("scalar_nudge_timescale", scalar_nudge_timescale),
        ("wind_nudge_timescale", wind_nudge_timescale),
    )
        τ > 0 || error("$name must be positive, got $τ s")
    end
    return SOCRATESForcing{FT, typeof(checked)}(
        checked,
        FT(1 / scalar_nudge_timescale),
        FT(1 / wind_nudge_timescale),
    )
end

Base.broadcastable(x::SOCRATESForcing) = tuple(x)

_t_seconds(t::Number) = Float64(t)
_t_seconds(t::CA.ClimaUtilities.TimeManager.ITime) = Float64(float(t))

"""
    sample_forcing(source, z, thermo)

`(; interpolants, surface)` for `source` on levels `z`: a per-level vector of time
interpolants for each of [`SSCF_FORCING_VARS`](@ref), and the surface series.
"""
function sample_forcing(c::SOCRATESCase, z::AbstractVector, thermo)
    interpolants = SSCF.get_column_forcing(
        c.flight_number,
        c.forcing_type,
        SSCF_FORCING_VARS;
        new_z = collect(Float64, z),
        thermodynamics_backend = thermo,
    )
    for name in SSCF_FORCING_VARS
        length(interpolants[name]) == length(z) || error(
            "SSCF returned $(length(interpolants[name])) level interpolants for \
             `$name` on a $(length(z))-level grid.",
        )
    end
    surface = SSCF.get_surface_forcing(
        c.flight_number,
        c.forcing_type;
        thermodynamics_backend = thermo,
    )
    return (; interpolants, surface)
end

sample_forcing(path::AbstractString, z::AbstractVector, thermo) =
    read_forcing_file(path, z, thermo)

# A ClimaCore broadcast is pointwise over the space and never exposes the level
# index needed to select `level_interpolants[k]`, so the values are written level
# by level. The guard requires one level per interpolant and one column: for a
# wider geometry `length(field)` exceeds the level count, and silently writing one
# column's profile into a flat index would be wrong.
function _fill_column!(field, level_interpolants, t::Float64)
    n = length(level_interpolants)
    (CC.Spaces.nlevels(axes(field)) == n && length(field) == n) || error(
        "SOCRATES forcing fills a single column of $n levels; got a field with \
         $(CC.Spaces.nlevels(axes(field))) levels and $(length(field)) points.",
    )
    FT = eltype(field)
    @inbounds for k in 1:n
        field[k] = FT(level_interpolants[k](t))
    end
    return nothing
end

"""
    ClimaAtmos.external_forcing_cache(Y, forcing::SOCRATESForcing, params, start_date)

Sample the forcing source onto the model levels and allocate the fields the
tendency writes.

Thermodynamic parameters come from `params`, the simulation's own set, so the
forcing and the model cannot disagree on the thermodynamics used to derive it.
"""
function CA.external_forcing_cache(Y, forcing::SOCRATESForcing, params, start_date)
    FT = CC.Spaces.undertype(axes(Y.c))
    thermo = CA.Parameters.thermodynamics_params(params)
    z = vec(Array(parent(CC.Fields.coordinate_field(Y.c).z)))
    (; interpolants, surface) = sample_forcing(forcing.source, z, thermo)
    return (;
        interpolants,
        surface,
        ᶠts = similar(CC.Fields.level(Y.f.u₃, CC.Utilities.half), FT),
        ᶜdTdt_hadv = similar(Y.c, FT),
        ᶜdqtdt_hadv = similar(Y.c, FT),
        ᶜT_nudge = similar(Y.c, FT),
        ᶜqt_nudge = similar(Y.c, FT),
        ᶜu_nudge = similar(Y.c, FT),
        ᶜv_nudge = similar(Y.c, FT),
        ᶜls_subsidence = similar(Y.c, FT),
    )
end

"""
    ClimaAtmos.external_forcing_tendency!(Yₜ, Y, p, t, forcing::SOCRATESForcing)

Apply the SOCRATES forcing at time `t`, seconds since the case start.

The interpolants are evaluated here rather than refreshed by a callback: they are
in-memory and exact at any `t`, so there is no staleness window. Horizontal
advection and scalar relaxation are summed into two scratch buffers and converted
once to `ρe_tot`/`ρq_tot` tendencies; wind relaxation and subsidence act on the
state directly.
"""
function CA.external_forcing_tendency!(Yₜ, Y, p, t, forcing::SOCRATESForcing)
    (;
        interpolants,
        ᶜdTdt_hadv,
        ᶜdqtdt_hadv,
        ᶜT_nudge,
        ᶜqt_nudge,
        ᶜu_nudge,
        ᶜv_nudge,
        ᶜls_subsidence,
    ) = p.external_forcing
    ts = _t_seconds(t)

    _fill_column!(ᶜdTdt_hadv, interpolants.dTdt_hadv, ts)
    _fill_column!(ᶜdqtdt_hadv, interpolants.dqtdt_hadv, ts)
    _fill_column!(ᶜT_nudge, interpolants.T_nudge, ts)
    _fill_column!(ᶜqt_nudge, interpolants.qt_nudge, ts)
    _fill_column!(ᶜu_nudge, interpolants.u_nudge, ts)
    _fill_column!(ᶜv_nudge, interpolants.v_nudge, ts)
    _fill_column!(ᶜls_subsidence, interpolants.subsidence, ts)

    CA.nudge_uv!(Yₜ, Y, p, ᶜu_nudge, ᶜv_nudge, forcing.inv_τ_wind)

    # `ᶜtemp_scalar`/`ᶜtemp_scalar_2` are the two slots ClimaAtmos's own forcing
    # tendencies use for this; the rest of the working fields are private cache
    # entries because scratch is shared across tendencies.
    ᶜdTdt = p.scratch.ᶜtemp_scalar
    ᶜdqtdt = p.scratch.ᶜtemp_scalar_2
    CA.nudge_Tq!(ᶜdTdt, ᶜdqtdt, Y, p, ᶜT_nudge, ᶜqt_nudge, forcing.inv_τ_scalar)
    @. ᶜdTdt += ᶜdTdt_hadv
    @. ᶜdqtdt += ᶜdqtdt_hadv
    CA.apply_Tq_forcing!(Yₜ, Y, p, ᶜdTdt, ᶜdqtdt)

    CA.apply_subsidence_forcing!(Yₜ, Y, p, ᶜls_subsidence)
    return nothing
end
