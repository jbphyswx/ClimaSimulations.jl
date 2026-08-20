"""
    model.jl

In-memory SSCF forcing, AtmosModel/Simulation assembly, diagnostic registration, and runner.
"""

using Dates: Dates
using ClimaAtmos: ClimaAtmos as CA
using ClimaComms: ClimaComms
using ClimaParams: ClimaParams as CP
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF

# --- In-Memory Column TimeVaryingInput ------------------------------------------------------ #

"""
    ColumnMemoryTimeVaryingInput(times, data, target_space; method)

A `TimeVaryingInput` over `data`, a `(n_levels, n_times)` array already on the target space's
levels, sampled at `times` (seconds since the simulation start date, sorted ascending).
"""
struct ColumnMemoryTimeVaryingInput{T <: AbstractVector, D <: AbstractMatrix, M} <:
       CA.ClimaUtilities.TimeVaryingInputs.AbstractTimeVaryingInput
    times::T
    data::D
    method::M

    function ColumnMemoryTimeVaryingInput(
        times::T,
        data::D,
        method::M,
    ) where {T <: AbstractVector, D <: AbstractMatrix, M}
        return new{T, D, M}(times, data, method)
    end
end

function ColumnMemoryTimeVaryingInput(
    times::AbstractVector,
    data::AbstractMatrix,
    target_space::CA.CC.Spaces.AbstractSpace;
    method = CA.ClimaUtilities.TimeVaryingInputs.LinearInterpolation(),
)
    issorted(times) || error("ColumnMemoryTimeVaryingInput times must be sorted ascending")
    length(times) >= 2 ||
        error("ColumnMemoryTimeVaryingInput needs at least 2 times, got $(length(times))")
    FT = CA.CC.Spaces.undertype(target_space)
    n_levels = _n_levels(target_space)
    size(data, 1) == n_levels || error(
        "ColumnMemoryTimeVaryingInput data has $(size(data, 1)) levels but the target space \
         has $n_levels. The data must already be sampled on the model's levels.",
    )
    size(data, 2) == length(times) || error(
        "ColumnMemoryTimeVaryingInput data has $(size(data, 2)) time samples but $(length(times)) \
         times were given",
    )
    array_type = ClimaComms.array_type(ClimaComms.device(target_space))
    t = collect(Float64, times)
    d = array_type(FT.(data))
    return ColumnMemoryTimeVaryingInput{typeof(t), typeof(d), typeof(method)}(t, d, method)
end

_n_levels(space) = size(parent(CA.CC.Fields.coordinate_field(space).z), 1)

Base.in(t::Number, itp::ColumnMemoryTimeVaryingInput) =
    first(itp.times) <= t <= last(itp.times)

_seconds(t::Number) = Float64(t)
_seconds(t::CA.ClimaUtilities.TimeManager.ITime) = Float64(float(t))

function CA.ClimaUtilities.TimeVaryingInputs.evaluate!(
    dest,
    itp::ColumnMemoryTimeVaryingInput,
    t,
    args...;
    kwargs...,
)
    ts = _seconds(t)
    i, w = _bracket(itp.times, ts, CA.ClimaUtilities.TimeVaryingInputs.extrapolation_bc(itp.method))
    dv = parent(dest)
    d = itp.data
    if w == 0
        @inbounds @views dv[:, 1] .= d[:, i]
    else
        @inbounds @views dv[:, 1] .= (1 - w) .* d[:, i] .+ w .* d[:, i + 1]
    end
    return nothing
end

function _bracket(times, ts, bc)
    if ts <= first(times)
        ts < first(times) && _check_in_range(bc, ts, times)
        return (firstindex(times), 0.0)
    elseif ts >= last(times)
        ts > last(times) && _check_in_range(bc, ts, times)
        return (lastindex(times), 0.0)
    end
    i = searchsortedlast(times, ts)
    t0, t1 = times[i], times[i + 1]
    return (i, (ts - t0) / (t1 - t0))
end

_check_in_range(::CA.ClimaUtilities.TimeVaryingInputs.Flat, ts, times) = nothing
_check_in_range(bc, ts, times) = error(
    "ColumnMemoryTimeVaryingInput evaluated at t = $ts s, outside its data range \
     [$(first(times)), $(last(times))] s.",
)

# --- In-Memory ColumnDataset Format --------------------------------------------------------- #

"""
    SocratesColumnFormat

Column-forcing format whose data lives in memory.
"""
struct SocratesColumnFormat <: CA.ColumnDatasets.AbstractColumnFormat end

CA.ColumnDatasets.format_name(::SocratesColumnFormat) = "SOCRATES (in-memory SSCF)"
CA.ColumnDatasets.format_variable_name(::SocratesColumnFormat, name::Symbol) = name
CA.ColumnDatasets.open_dataset(f, ::SocratesColumnFormat, path, options) = f(options)
CA.ColumnDatasets.height_profile(::SocratesColumnFormat, ds, options) = ds.z
CA.ColumnDatasets.dates(::SocratesColumnFormat, ds) =
    ds.start_date .+ Dates.Millisecond.(round.(Int, 1000 .* ds.times))
CA.ColumnDatasets.has_variable(::SocratesColumnFormat, ds, name::Symbol) =
    haskey(ds.column, name) || haskey(ds.surface, name)
CA.ColumnDatasets.read_profile(::SocratesColumnFormat, ds, name::Symbol, time_index) =
    ds.column[name][:, time_index]
CA.ColumnDatasets.read_series(::SocratesColumnFormat, ds, name::Symbol) = ds.surface[name]
CA.ColumnDatasets.site_location(::SocratesColumnFormat, ds) =
    (; latitude = ds.lat, longitude = ds.lon)

function CA.ColumnDatasets.column_timevaryinginputs(
    cd::CA.ColumnDatasets.ColumnDataset{SocratesColumnFormat},
    names,
    target_space,
    start_date;
    method = CA.ColumnDatasets.time_interpolation_method(cd.format),
)
    names = Tuple(names)
    inputs = CA.ColumnDatasets.open_dataset(cd) do ds
        start_date == ds.start_date || error(
            "SOCRATES forcing was sampled with start_date $(ds.start_date) but the simulation \
             uses $start_date.",
        )
        map(names) do name
            haskey(ds.column, name) || error(
                "SOCRATES forcing has no column variable `$name`; it provides \
                 $(sort(collect(keys(ds.column)))).",
            )
            ColumnMemoryTimeVaryingInput(ds.times, ds.column[name], target_space; method)
        end
    end
    return NamedTuple{names}(inputs)
end

# --- Forcing Sampling & Memoization -------------------------------------------------------- #

const SSCF_TO_CANONICAL = (
    dTdt_hadv = :tntha,
    dqtdt_hadv = :tnhusha,
    T_nudge = :ta,
    qt_nudge = :hus,
    u_nudge = :ua,
    v_nudge = :va,
    subsidence = :wa,
)

const SSCF_COLUMN_KEYS = Tuple(keys(SSCF_TO_CANONICAL))
const DEFAULT_FORCING_DT = 300.0
const _FORCING_CACHE = Dict{Any, Any}()

_forcing_cache_key(FT, c, z, dt_sec, start_date) =
    (FT, c.flight_number, forcing_label(c), hash(z), dt_sec, start_date)

"""
    socrates_forcing_arrays(FT, case; z, dt_sec, start_date, refresh)

Sample SSCF forcing for `case` onto levels `z` and time axis `dt_sec`. Memoized.
"""
function socrates_forcing_arrays(
    ::Type{FT},
    c::SocratesCase;
    z::AbstractVector,
    dt_sec::Real = DEFAULT_FORCING_DT,
    start_date::Dates.DateTime = simulation_start_date(c),
    refresh::Bool = false,
) where {FT <: AbstractFloat}
    key = _forcing_cache_key(FT, c, z, dt_sec, start_date)
    if !refresh && haskey(_FORCING_CACHE, key)
        return _FORCING_CACHE[key]
    end
    arrays = _build_forcing_arrays(FT, c; z, dt_sec, start_date)
    _FORCING_CACHE[key] = arrays
    return arrays
end

forcing_cache_size() = length(_FORCING_CACHE)
empty_forcing_cache!() = (empty!(_FORCING_CACHE); nothing)

function _build_forcing_arrays(
    ::Type{FT},
    c::SocratesCase;
    z::AbstractVector,
    dt_sec::Real,
    start_date::Dates.DateTime,
) where {FT <: AbstractFloat}
    validate(c)
    ft = c.forcing_type
    flight = c.flight_number
    z_model = collect(FT, z)
    thermo = _thermodynamics_backend(FT)

    times = collect(FT, 0:FT(dt_sec):(FT(t_end(c)) + FT(dt_sec)))

    column = SSCF.get_column_forcing(
        flight,
        ft,
        SSCF_COLUMN_KEYS;
        new_z = z_model,
        thermodynamics_backend = thermo,
    )
    surface = SSCF.get_surface_forcing(flight, ft; thermodynamics_backend = thermo)

    sampled = Dict{Symbol, Matrix{FT}}()
    for (sscf_name, canonical) in pairs(SSCF_TO_CANONICAL)
        sampled[canonical] = _sample_levels(FT, column[sscf_name], z_model, times)
    end
    sampled[:tntva] = zeros(FT, length(z_model), length(times))
    sampled[:tnhusva] = zeros(FT, length(z_model), length(times))
    sampled[:rho] = _initial_density(FT, c, z_model, length(times))

    lat, lon = _site_location(FT, c)
    coszen, rsdt = _insolation_series(FT, c, times, lat, lon)
    surface_series = Dict{Symbol, Vector{FT}}(
        :ts => FT[FT(surface.Tsfc(t)) for t in times],
        :coszen => coszen,
        :rsdt => rsdt,
    )

    return (;
        z = z_model,
        times,
        start_date,
        lat,
        lon,
        column = sampled,
        surface = surface_series,
    )
end

function _sample_levels(
    ::Type{FT},
    level_interpolants,
    z_model::AbstractVector,
    times::AbstractVector,
) where {FT}
    length(level_interpolants) == length(z_model) || error(
        "SSCF returned $(length(level_interpolants)) level interpolants for a $(length(z_model))-level grid.",
    )
    out = Matrix{FT}(undef, length(z_model), length(times))
    for (j, t) in enumerate(times), k in eachindex(z_model)
        out[k, j] = FT(level_interpolants[k](t))
    end
    return out
end

function _initial_density(
    ::Type{FT},
    c::SocratesCase,
    z_model::AbstractVector,
    n_times::Int,
) where {FT}
    les = SSCF.open_atlas_les_output(c.flight_number, c.forcing_type)
    z_les = FT.(vec(Array(les.data["z"])))
    ρ_les = FT.(Array(les.data["RHO"])[:, 1])
    ρ = _linear_interp(z_model, z_les, ρ_les)
    return repeat(reshape(ρ, :, 1), 1, n_times)
end

function _linear_interp(x_out::AbstractVector, xs::AbstractVector, ys::AbstractVector)
    issorted(xs) || error("_linear_interp requires ascending source coordinates")
    FT = eltype(ys)
    out = similar(x_out, FT)
    for (i, x) in enumerate(x_out)
        if x <= first(xs)
            out[i] = first(ys)
        elseif x >= last(xs)
            out[i] = last(ys)
        else
            j = searchsortedlast(xs, x)
            w = (x - xs[j]) / (xs[j + 1] - xs[j])
            out[i] = (1 - w) * ys[j] + w * ys[j + 1]
        end
    end
    return out
end

function _site_location(::Type{FT}, c::SocratesCase) where {FT}
    inp = SSCF.open_atlas_les_input(c.flight_number, c.forcing_type)
    lat = FT(Array(inp.data["lat"])[1])
    lon = FT(Array(inp.data["lon"])[1])
    return (FT(lat), FT(mod(lon + 180, 360) - 180))
end

function _insolation_series(
    ::Type{FT},
    c::SocratesCase,
    times::AbstractVector,
    lat,
    lon,
) where {FT}
    params = CA.Insolation.Parameters.InsolationParameters(FT)
    reference = les_start_datetime(c) + Dates.Hour(12)
    F, _, μ, _ = CA.Insolation.insolation(reference, FT(lat), FT(lon), params)
    return fill(FT(μ), length(times)), fill(FT(F), length(times))
end

_thermodynamics_backend(::Type{FT}) where {FT} =
    CA.Thermodynamics.Parameters.ThermodynamicsParameters(CP.create_toml_dict(FT))

const SCALAR_NUDGE_TIMESCALE = 20.0 * 60.0
const OBS_WIND_NUDGE_TIMESCALE = 20.0 * 60.0
const ERA5_WIND_NUDGE_TIMESCALE = 60.0 * 60.0

wind_nudge_timescale(c::SocratesCase) =
    c.forcing_type isa SSCF.ObsForcing ? OBS_WIND_NUDGE_TIMESCALE :
    ERA5_WIND_NUDGE_TIMESCALE

default_socrates_forcing_terms(c::SocratesCase) = (
    CA.HorizontalAdvection(),
    CA.Nudging(:ta, :hus; timescale = SCALAR_NUDGE_TIMESCALE),
    CA.Nudging(:ua, :va; timescale = wind_nudge_timescale(c)),
    CA.Subsidence(),
)

"""
    socrates_forcing(FT, case; z, dt_sec, start_date, forcing_terms, refresh)

The `ExternalDrivenTVForcing` for `case`, backed by the in-memory dataset.
"""
function socrates_forcing(
    ::Type{FT},
    c::SocratesCase;
    z::AbstractVector,
    dt_sec::Real = DEFAULT_FORCING_DT,
    start_date::Dates.DateTime = simulation_start_date(c),
    forcing_terms = default_socrates_forcing_terms(c),
    refresh::Bool = false,
) where {FT <: AbstractFloat}
    arrays = socrates_forcing_arrays(FT, c; z, dt_sec, start_date, refresh)
    dataset = CA.ColumnDatasets.ColumnDataset(
        "sscf:$(case_name(c))";
        format = SocratesColumnFormat(),
        arrays...,
    )
    return CA.ExternalDrivenTVForcing(dataset; forcing = forcing_terms)
end

# --- Scaled Terminal Velocity Mode ---------------------------------------------------------- #

"""
    ScaledTerminalVelocity{FT}(liquid, ice, rain, snow)

Diagnostic size-dependent terminal velocity scaled by one factor per species.
"""
struct ScaledTerminalVelocity{FT} <: CA.AbstractTerminalVelocityMode
    liquid::FT
    ice::FT
    rain::FT
    snow::FT
end

function terminal_velocity_scaling(::Type{FT}, toml_dict) where {FT <: AbstractFloat}
    factors = CP.get_parameter_values(
        toml_dict,
        TERMINAL_VELOCITY_SCALING_PARAMS,
        "SOCRATES",
    )
    return ScaledTerminalVelocity{FT}(
        factors.liquid,
        factors.ice,
        factors.rain,
        factors.snow,
    )
end

for (field, tracer) in
    ((:liquid, :q_lcl), (:ice, :q_icl), (:rain, :q_rai), (:snow, :q_sno))
    @eval CA.terminal_velocity(
        model::CA.NonEquilibriumMicrophysics1M,
        mode::ScaledTerminalVelocity,
        name::CA.CC.MatrixFields.FieldName{($(QuoteNode(tracer)),)},
        cmc,
        cmp,
        ρ,
        q,
    ) =
        mode.$field * CA.terminal_velocity(
            model,
            CA.DiagnosticTerminalVelocity(),
            name,
            cmc,
            cmp,
            ρ,
            q,
        )
end

CA.gs_terminal_velocity(
    model::CA.NonEquilibriumMicrophysics1M,
    ::ScaledTerminalVelocity,
    name,
    ρwχ,
    ρχ,
) = CA.gs_terminal_velocity(
    model,
    CA.DiagnosticTerminalVelocity(),
    name,
    ρwχ,
    ρχ,
)

# --- Diagnostics Registration -------------------------------------------------------------- #

const PROFILE_VARS = ("clw", "cli", "husra", "hussn")
const PATH_VARS = ("lwp", "iwp", "rwp", "swp")
const DEFAULT_DIAGNOSTIC_VARS = (PROFILE_VARS..., PATH_VARS...)
const DEFAULT_DIAGNOSTIC_PERIOD = "10mins"

const TRANSPORT_SPECIES = (
    ("q_lcl", :ρq_lcl, :ᶜwₗ),
    ("q_icl", :ρq_icl, :ᶜwᵢ),
    ("q_rai", :ρq_rai, :ᶜwᵣ),
    ("q_sno", :ρq_sno, :ᶜwₛ),
)

function sedimentation_tendency(state, cache, ρq_name::Symbol, w_name::Symbol)
    ᶜρ = state.c.ρ
    FT = eltype(ᶜρ)
    ᶜJ = CA.CC.Fields.local_geometry_field(axes(state.c)).J
    ᶠJ = CA.CC.Fields.local_geometry_field(axes(state.f)).J
    ᶜw = getproperty(cache.precomputed, w_name)
    ᶜρq = getproperty(state.c, ρq_name)
    ᶠinterp = CA.CC.Operators.InterpolateC2F()
    ᶠright_bias = CA.CC.Operators.RightBiasedC2F()
    ᶜprecipdivᵥ = CA.CC.Operators.DivergenceF2C(
        top = CA.CC.Operators.SetValue(CA.CC.Geometry.Contravariant3Vector(zero(FT))),
    )
    return @. CA.LazyBroadcast.lazy(
        -(ᶜprecipdivᵥ(
            ᶠinterp(ᶜρ * ᶜJ) / ᶠJ *
            ᶠright_bias(CA.CC.Geometry.WVector(-(ᶜw)) * CA.specific(ᶜρq, ᶜρ)),
        )) / ᶜρ,
    )
end

function vertical_advection_tendency(state, cache, ρq_name::Symbol)
    ᶜρ = state.c.ρ
    ᶠu³ = cache.precomputed.ᶠu³
    ᶜρq = getproperty(state.c, ρq_name)
    ᶜχ = @. CA.LazyBroadcast.lazy(CA.specific(ᶜρq, ᶜρ))
    vtt = CA.vertical_transport(
        ᶜρ,
        ᶠu³,
        ᶜχ,
        cache.dt,
        cache.atmos.numerics.tracer_upwinding,
    )
    return @. CA.LazyBroadcast.lazy(vtt / ᶜρ)
end

function vertical_diffusion_tendency(state, cache, ρq_name::Symbol)
    ᶜρ = state.c.ρ
    FT = eltype(ᶜρ)
    (; ᶠK_h, ᶠK_entr) = cache.precomputed
    α = CA.CAP.α_vert_diff_tracer(cache.params)
    ᶜρq = getproperty(state.c, ρq_name)
    ᶜχ = @. CA.LazyBroadcast.lazy(CA.specific(ᶜρq, ᶜρ))
    ᶠgradᵥ = CA.CC.Operators.GradientC2F()
    ᶜdivᵥ_ρq = CA.CC.Operators.DivergenceF2C(
        top = CA.CC.Operators.SetValue(CA.CC.Geometry.Covariant3Vector(zero(FT))),
        bottom = CA.CC.Operators.SetValue(CA.CC.Geometry.Covariant3Vector(zero(FT))),
    )
    ᶠinterp = CA.CC.Operators.InterpolateC2F()
    ᶠρK = @. CA.LazyBroadcast.lazy(ᶠinterp(ᶜρ) * (α * ᶠK_h + ᶠK_entr))
    return @. CA.LazyBroadcast.lazy(ᶜdivᵥ_ρq(ᶠρK * ᶠgradᵥ(ᶜχ)) / ᶜρ)
end

function register_transport_diagnostics!()
    for (species, ρq_name, w_name) in TRANSPORT_SPECIES
        for (prefix, description, compute) in (
            ("sed", "sedimentation", (state, cache, _) -> sedimentation_tendency(state, cache, ρq_name, w_name)),
            ("adv", "resolved vertical advection", (state, cache, _) -> vertical_advection_tendency(state, cache, ρq_name)),
            ("dif", "SGS vertical diffusion", (state, cache, _) -> vertical_diffusion_tendency(state, cache, ρq_name)),
        )
            name = "$(prefix)_$(species)"
            delete!(CA.Diagnostics.ALL_DIAGNOSTICS, name)
            CA.Diagnostics.add_diagnostic_variable!(;
                short_name = name,
                long_name = "Tendency of $species from $description",
                units = "kg kg^-1 s^-1",
                comments = "Computed with the model's own vertical operators on the model state.",
                compute,
            )
        end
    end
    return nothing
end

function socrates_diagnostics(
    short_names = DEFAULT_DIAGNOSTIC_VARS;
    period::AbstractString = DEFAULT_DIAGNOSTIC_PERIOD,
    reduction::AbstractString = "average",
    n_levels::Int,
    output_at_levels::Bool = true,
)
    register_transport_diagnostics!()
    return CA.DiagnosticsConfig(;
        default = false,
        additional = ((; short_name = collect(String.(short_names)), period, reduction),),
        interpolation_num_points = (2, 2, n_levels),
        output_at_levels,
    )
end

# --- Setup & Model Construction ------------------------------------------------------------ #

const DEFAULT_Z0 = 1.0e-4

struct SocratesSetup{S}
    inner::S
end

function socrates_setup(
    ::Type{FT},
    c::SocratesCase;
    z::AbstractVector,
    dt_sec::Real = DEFAULT_FORCING_DT,
    start_date::Dates.DateTime = simulation_start_date(c),
    forcing_terms = default_socrates_forcing_terms(c),
    z0::Real = DEFAULT_Z0,
    refresh::Bool = false,
) where {FT <: AbstractFloat}
    forcing = socrates_forcing(FT, c; z, dt_sec, start_date, forcing_terms, refresh)
    return SocratesSetup(
        CA.Setups.ForcingFromFile(
            forcing.dataset,
            Dates.format(start_date, "yyyymmdd");
            forcing,
            flux_scheme = CA.SurfaceConditions.MoninObukhov(; z0 = FT(z0)),
        ),
    )
end

function CA.Setups.center_initial_condition(
    setup::SocratesSetup,
    local_geometry,
    params,
)
    profiles = setup.inner.profiles
    (; z) = local_geometry.coordinates
    FT = typeof(z)
    T = FT(profiles.T(z))
    ρ = FT(profiles.ρ(z))
    q_tot = FT(profiles.q_tot(z))
    thermo = CA.Parameters.thermodynamics_params(params)
    q_sat_liq = CA.TD.q_vap_saturation(thermo, T, ρ, CA.TD.Liquid())
    return CA.Setups.physical_state(;
        T,
        ρ,
        q_tot,
        q_liq = max(zero(FT), q_tot - FT(q_sat_liq)),
        q_ice = zero(FT),
        u = FT(profiles.u(z)),
        v = FT(profiles.v(z)),
        tke = zero(FT),
    )
end

CA.Setups.surface_condition(setup::SocratesSetup, params) =
    CA.Setups.surface_condition(setup.inner, params)

CA.Setups.external_forcing(setup::SocratesSetup, ::Type{FT}) where {FT} =
    CA.Setups.external_forcing(setup.inner, FT)

CA.Setups.insolation_model(setup::SocratesSetup) =
    CA.Setups.insolation_model(setup.inner)

CA.Setups.surface_temperature_model(setup::SocratesSetup) =
    CA.Setups.surface_temperature_model(setup.inner)

function socrates_model(
    ::Type{FT},
    params;
    external_forcing,
    area_fraction::Real = 1.0e-3,
    z0::Real = DEFAULT_Z0,
    surface_albedo::Real = 0.07,
    terminal_velocity_mode = ScaledTerminalVelocity{FT}(one(FT), one(FT), one(FT), one(FT)),
    edmfx = CA.EDMFXModel(;
        entr_model = CA.PiGroupsEntrainment(),
        detr_model = CA.BuoyancyVelocityDetrainment(),
        sgs_mass_flux = true,
        sgs_diffusive_flux = true,
        nh_pressure = true,
        vertical_diffusion = true,
        filter = true,
        scale_blending_method = CA.SmoothMinimumBlending(),
    ),
    turbconv = CA.AtmosTurbconv(;
        edmfx_model = edmfx,
        turbconv_model = CA.PrognosticEDMFX(;
            n_updrafts = 1,
            prognostic_tke = true,
            area_fraction = FT(area_fraction),
        ),
    ),
    water = CA.AtmosWater(;
        microphysics_model = CA.NonEquilibriumMicrophysics1M(;
            n_substeps = 3,
            n_substeps_quad = 2,
        ),
        cloud_model = CA.QuadratureCloud(),
        microphysics_tendency_timestepping = CA.Implicit(),
        tracer_nonnegativity_method = nothing,
        sgs_quadrature = CA.SGSQuadrature(
            FT;
            quadrature_order = 3,
            distribution = CA.GaussianSGS(),
            T_min = FT(CA.Parameters.T_min_sgs(params)),
            q_max = FT(CA.Parameters.q_max_sgs(params)),
        ),
        terminal_velocity_mode,
    ),
    radiation = CA.AtmosRadiation(;
        radiation_mode = CA.RRTMGPInterface.AllSkyRadiationWithClearSkyDiagnostics(),
        insolation = CA.ExternalTVInsolation(),
    ),
    surface = CA.AtmosSurface(;
        flux_scheme = CA.SurfaceConditions.MoninObukhov(; z0 = FT(z0)),
        temperature = CA.SurfaceConditions.ExternalTemperature(),
        surface_albedo = CA.ConstantAlbedo{FT}(; α = FT(surface_albedo)),
    ),
    numerics = CA.AtmosNumerics(;
        diff_mode = CA.Implicit(),
        hyperdiff = nothing,
        edmfx_sgsflux_upwinding = :first_order,
    ),
    sponge = CA.AtmosSponge(),
) where {FT <: AbstractFloat}
    return CA.AtmosModel(;
        water,
        scm_setup = CA.SCMSetup(; external_forcing),
        radiation,
        turbconv,
        surface,
        numerics,
        sponge,
    )
end

function socrates_ode_config(
    ::Type{FT};
    ode_algo::AbstractString = "ARS222",
    update_jacobian_every::AbstractString = "stage",
    max_newton_iters::Int = 1,
    use_krylov_method = false,
    use_dynamic_krylov_rtol = false,
    eisenstat_walker_forcing_alpha = 2.0,
    krylov_rtol = 0.1,
    use_newton_rtol = false,
    newton_rtol = 1.0e-5,
    jvp_step_adjustment = 1.0,
) where {FT <: AbstractFloat}
    return CA.ode_configuration(
        FT,
        ode_algo,
        update_jacobian_every,
        max_newton_iters,
        use_krylov_method,
        use_dynamic_krylov_rtol,
        eisenstat_walker_forcing_alpha,
        krylov_rtol,
        use_newton_rtol,
        newton_rtol,
        jvp_step_adjustment,
    )
end

"""
    socrates_simulation(FT, case; kwargs...)

Assemble an `AtmosSimulation{FT}` for `case`.
"""
function socrates_simulation(
    ::Type{FT},
    c::SocratesCase;
    params = nothing,
    grid = socrates_grid(FT, c),
    t_end::Real = t_end(c),
    t_start::Real = 0,
    dt::Real = 10,
    output_dir::AbstractString,
    diagnostics = nothing,
    forcing_dt::Real = DEFAULT_FORCING_DT,
    forcing_terms = default_socrates_forcing_terms(c),
    ode_config = socrates_ode_config(FT),
    jacobian = CA.ManualSparseJacobian(; approximate_solve_iters = 2),
    checkpoint_frequency = Inf,
    area_fraction::Real = 1.0e-3,
    z0::Real = DEFAULT_Z0,
    model_kwargs = (;),
    job_id::AbstractString = case_name(c),
    output_dir_style::AbstractString = "activelink",
    verbose::Bool = true,
    refresh_forcing::Bool = false,
) where {FT <: AbstractFloat}
    validate(c)
    start_date = simulation_start_date(c)
    z = socrates_z(grid)
    setup = socrates_setup(
        FT,
        c;
        z,
        dt_sec = forcing_dt,
        start_date,
        forcing_terms,
        z0,
        refresh = refresh_forcing,
    )
    toml_dict = socrates_toml_dict(FT, c; params)
    atmos_params = socrates_params(toml_dict)
    return CA.AtmosSimulation{FT}(;
        model = socrates_model(
            FT,
            atmos_params;
            external_forcing = CA.Setups.external_forcing(setup, FT),
            area_fraction,
            z0,
            terminal_velocity_mode = terminal_velocity_scaling(FT, toml_dict),
            model_kwargs...,
        ),
        grid,
        setup,
        params = atmos_params,
        dt,
        t_start,
        t_end,
        start_date,
        ode_config,
        jacobian,
        diagnostics = something(diagnostics, socrates_diagnostics(; n_levels = length(z))),
        job_id,
        output_dir,
        output_dir_style,
        checkpoint_frequency,
        verbose,
    )
end

# Alias
const simulation = socrates_simulation

# --- Execution Infrastructure -------------------------------------------------------------- #

abstract type AbstractExecutor end

struct SerialExecutor <: AbstractExecutor end

function run_tasks(f, tasks, ::SerialExecutor = SerialExecutor())
    results = Vector{Any}(nothing, length(tasks))
    for (i, task) in enumerate(tasks)
        try
            results[i] = f(task)
        catch e
            @error "Task failed" task exception = (e, catch_backtrace())
        end
    end
    return results
end

"""
    run_case(case; FT, output_dir, kwargs...) -> output_dir

Build and solve one SOCRATES case, returning the directory its diagnostics were written to.
"""
function run_case(
    c::SocratesCase;
    FT::Type{<:AbstractFloat} = Float64,
    output_dir::AbstractString,
    kwargs...,
)
    sim = socrates_simulation(FT, c; output_dir, kwargs...)
    result = CA.solve_atmos!(sim)
    result.ret_code === :success || error(
        "$(case_name(c)) did not run to completion: ClimaAtmos returned `$(result.ret_code)`.",
    )
    return sim.output_dir
end

"""
    run_cases(cases; FT, output_dir, executor, kwargs...) -> Vector

Run multiple cases.
"""
function run_cases(
    cases::AbstractVector{<:SocratesCase};
    FT::Type{<:AbstractFloat} = Float64,
    output_dir::AbstractString,
    executor::AbstractExecutor = SerialExecutor(),
    kwargs...,
)
    foreach(validate, cases)
    run_one(c) = run_case(c; FT, output_dir = joinpath(output_dir, case_name(c)), kwargs...)
    return run_tasks(run_one, cases, executor)
end
