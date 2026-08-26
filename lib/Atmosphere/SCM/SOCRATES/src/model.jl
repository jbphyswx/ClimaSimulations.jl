"""
    model.jl

Model components SOCRATES owns — scaled sedimentation, prescribed surface
temperature and insolation, the initial condition — and the assembly of an
`AtmosSimulation` from them.

Every default here is overridable: `socrates_model` and `socrates_simulation`
forward `kwargs...` to `AtmosModel` and `AtmosSimulation`, so any field of either
can be replaced without editing this file.
"""

# --- Sedimentation velocity scaling ----------------------------------------- #

"""
    ScaledTerminalVelocity{FT}(scale)

Diagnostic terminal velocity of one species, scaled by `scale`. A scale of 1
reproduces `ClimaAtmos.DiagnosticTerminalVelocity` exactly.
"""
struct ScaledTerminalVelocity{FT} <: CA.AbstractTerminalVelocityMode
    scale::FT
end

CA.terminal_velocity(
    model::CA.NonEquilibriumMicrophysics1M,
    mode::ScaledTerminalVelocity,
    name,
    args...,
) =
    mode.scale *
    CA.terminal_velocity(model, CA.DiagnosticTerminalVelocity(), name, args...)

# `scale` is deliberately not reapplied here: the diagnostic method forms
# `ρwχ / ρχ` from subdomain fluxes that were already scaled above.
CA.gs_terminal_velocity(
    model::CA.NonEquilibriumMicrophysics1M,
    mode::ScaledTerminalVelocity,
    name,
    args...,
) = CA.gs_terminal_velocity(model, CA.DiagnosticTerminalVelocity(), name, args...)

"""Parameter name to species for the sedimentation scaling factors."""
const TERMINAL_VELOCITY_SCALING_PARAMS = (;
    :cloud_liquid_terminal_velocity_scaling_factor => :liquid,
    :cloud_ice_terminal_velocity_scaling_factor => :ice,
    :rain_terminal_velocity_scaling_factor => :rain,
    :snow_terminal_velocity_scaling_factor => :snow,
)

"""
    terminal_velocity_modes(FT, toml_dict)

The four `ScaledTerminalVelocity` modes, read from `toml_dict`.
"""
function terminal_velocity_modes(::Type{FT}, toml_dict) where {FT}
    scaling = CP.get_parameter_values(
        toml_dict,
        TERMINAL_VELOCITY_SCALING_PARAMS,
        "SOCRATES",
    )
    return (;
        terminal_velocity_liquid = ScaledTerminalVelocity{FT}(FT(scaling.liquid)),
        terminal_velocity_ice = ScaledTerminalVelocity{FT}(FT(scaling.ice)),
        terminal_velocity_rain = ScaledTerminalVelocity{FT}(FT(scaling.rain)),
        terminal_velocity_snow = ScaledTerminalVelocity{FT}(FT(scaling.snow)),
    )
end

# --- Prescribed surface temperature ----------------------------------------- #

"""
    SOCRATESSurfaceTemperature()

Sea surface temperature read from the case's SSCF surface forcing.

Requires a [`SOCRATESForcing`](@ref), whose cache carries the surface splines and
the field written here.
"""
struct SOCRATESSurfaceTemperature <: CA.SurfaceConditions.SurfaceTemperature end

function CA.SurfaceConditions.surface_temperature(
    ::SOCRATESSurfaceTemperature,
    Y,
    p,
    t_time,
)
    (; surface, ᶠts) = p.external_forcing
    fill!(parent(ᶠts), eltype(ᶠts)(surface.Tsfc(_t_seconds(t_time))))
    return CC.Fields.field_values(ᶠts)
end

# --- Prescribed insolation -------------------------------------------------- #

"""
    SOCRATESInsolation{FT}(cos_zenith, toa_flux)

Insolation held at a fixed zenith angle, following Atlas et al. (2020) section
4.2: "The solar zenith angle is held constant at the reference time of the case."

`toa_flux` is the solar flux perpendicular to the beam, so RRTMGP's downwelling
shortwave `toa_flux * cos_zenith` is the insolation on a horizontal surface.
"""
struct SOCRATESInsolation{FT} <: CA.AbstractInsolation
    cos_zenith::FT
    toa_flux::FT
end

"""
    SOCRATESInsolation(FT, flight_number, forcing_type; reference_time, latitude, longitude)

Evaluate the orbital insolation for a case. By default at the Atlas reference
time — 12 h after the LES start — and at the site of the flight.
"""
function SOCRATESInsolation(
    ::Type{FT},
    flight_number::Integer,
    forcing_type::SSCF.AbstractForcingType;
    reference_time::Dates.DateTime = SSCF.get_socrates_initial_time(
        Int(flight_number),
    ) + Dates.Hour(12),
    site = site_location(FT, flight_number, forcing_type),
    insolation_params = Insolation.Parameters.InsolationParameters(FT),
) where {FT <: AbstractFloat}
    latitude, longitude = site
    (; S, μ) =
        Insolation.insolation(reference_time, latitude, longitude, insolation_params)
    μ > 0 || error(
        "The reference time $reference_time for flight $flight_number is at night \
         (cos(zenith) = $μ). RRTMGP needs a positive zenith cosine; pass a \
         polar-night insolation model instead.",
    )
    return SOCRATESInsolation{FT}(FT(μ), FT(S))
end

function CA.set_insolation_variables!(Y, p, t, insolation::SOCRATESInsolation)
    (; rrtmgp_solver) = p.radiation
    CA.RRTMGP.cos_zenith(rrtmgp_solver) .= insolation.cos_zenith
    CA.RRTMGP.toa_flux(rrtmgp_solver) .= insolation.toa_flux
    return nothing
end

"""
    site_location(FT, flight_number, forcing_type)

`(latitude, longitude)` of the flight in degrees, longitude in `[-180, 180)`.
"""
function site_location(
    ::Type{FT},
    flight_number::Integer,
    forcing_type::SSCF.AbstractForcingType,
) where {FT}
    inp = SSCF.open_atlas_les_input(Int(flight_number), forcing_type)
    lat = FT(Array(inp.data["lat"])[1])
    lon = FT(Array(inp.data["lon"])[1])
    return (lat, FT(mod(lon + 180, 360) - 180))
end

# --- Initial condition ------------------------------------------------------ #

"""
    SOCRATESSetup(FT, flight_number, forcing_type, z; thermo_params, density)

Initial state of a SOCRATES case: the `t = 0` sounding profiles on levels `z`,
with the density the Atlas LES started from.
"""
struct SOCRATESSetup{P <: CA.Setups.ColumnProfiles}
    profiles::P
end

function SOCRATESSetup(
    ::Type{FT},
    flight_number::Integer,
    forcing_type::SSCF.AbstractForcingType,
    z::AbstractVector;
    thermo_params,
    density = initial_density(FT, flight_number, forcing_type, z),
) where {FT <: AbstractFloat}
    ic = SSCF.get_column_forcing(
        Int(flight_number),
        forcing_type,
        (:T_nudge, :qt_nudge, :u_nudge, :v_nudge);
        new_z = collect(FT, z),
        thermodynamics_backend = thermo_params,
        initial_condition = true,
    )
    return SOCRATESSetup(
        CA.Setups.ColumnProfiles(
            collect(FT, z),
            collect(FT, ic.T_nudge),
            collect(FT, ic.u_nudge),
            collect(FT, ic.v_nudge),
            collect(FT, ic.qt_nudge),
            collect(FT, density),
        ),
    )
end

"""
    initial_density(FT, flight_number, forcing_type, z)

The density [kg/m³] the Atlas LES started from, on levels `z`: its `RHO` at the
first output time, interpolated in height.
"""
function initial_density(
    ::Type{FT},
    flight_number::Integer,
    forcing_type::SSCF.AbstractForcingType,
    z::AbstractVector,
) where {FT}
    les = SSCF.open_atlas_les_output(Int(flight_number), forcing_type)
    z_les = FT.(vec(Array(les.data["z"])))
    ρ_les = FT.(Array(les.data["RHO"])[:, 1])
    interp = CA.Setups.ColumnProfiles(z_les, ρ_les, ρ_les, ρ_les, ρ_les, ρ_les).T
    return FT[interp(FT(zk)) for zk in z]
end

"""
    ClimaAtmos.Setups.center_initial_condition(setup::SOCRATESSetup, local_geometry, params)

The sounding state at one level, with any supersaturation placed in cloud liquid.

Without that split the column starts supersaturated and the microphysics deposits
the condensate over the first steps, releasing latent heat the observed `T`
already contains.
"""
function CA.Setups.center_initial_condition(
    setup::SOCRATESSetup,
    local_geometry,
    params,
)
    (; profiles) = setup
    (; z) = local_geometry.coordinates
    FT = typeof(z)
    T = FT(profiles.T(z))
    ρ = FT(profiles.ρ(z))
    q_tot = FT(profiles.q_tot(z))
    thermo = CA.Parameters.thermodynamics_params(params)
    q_sat_liq = TD.q_vap_saturation(thermo, T, ρ, TD.Liquid())
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

# --- Diagnostics ------------------------------------------------------------ #

"""The profile and water-path variables the Atlas LES comparison scores."""
const SCORED_VARS = ("clw", "cli", "husra", "hussn", "lwp", "iwp", "rwp", "swp")

"""
    socrates_diagnostics(short_names; period_seconds, reduction, n_levels, kwargs...)

A `DiagnosticsConfig` writing `short_names` on the model's own levels.
`kwargs...` are forwarded to `DiagnosticsConfig`.
"""
socrates_diagnostics(
    short_names = SCORED_VARS;
    period_seconds::Real = 600,
    reduction::AbstractString = "average",
    n_levels::Int,
    kwargs...,
) = CA.DiagnosticsConfig(;
    default = false,
    additional = ((;
        short_name = collect(String.(short_names)),
        period = "$(period_seconds)secs",
        reduction,
    ),),
    interpolation_num_points = (2, 2, n_levels),
    output_at_levels = true,
    kwargs...,
)

# --- Model ------------------------------------------------------------------ #

"""
    socrates_model(FT, params; external_forcing, insolation, kwargs...)

The `AtmosModel` for a SOCRATES case: prognostic EDMFX with one updraft,
non-equilibrium 1-moment microphysics with a quadrature cloud, interactive
RRTMGP radiation, and Monin-Obukhov surface fluxes over a prescribed SST.

`kwargs...` are forwarded to `AtmosModel` and override any default below, so any
component can be replaced without editing this function.
"""
function socrates_model(
    ::Type{FT},
    params;
    external_forcing,
    insolation,
    # water
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
    terminal_velocity_liquid = ScaledTerminalVelocity{FT}(one(FT)),
    terminal_velocity_ice = ScaledTerminalVelocity{FT}(one(FT)),
    terminal_velocity_rain = ScaledTerminalVelocity{FT}(one(FT)),
    terminal_velocity_snow = ScaledTerminalVelocity{FT}(one(FT)),
    # radiation
    radiation_mode = CA.RRTMGPInterface.AllSkyRadiationWithClearSkyDiagnostics(),
    # turbulence-convection
    edmfx_model = CA.EDMFXModel(;
        entr_model = CA.PiGroupsEntrainment(),
        detr_model = CA.BuoyancyVelocityDetrainment(),
        sgs_mass_flux = true,
        sgs_diffusive_flux = true,
        nh_pressure = true,
        vertical_diffusion = true,
        filter = true,
        scale_blending_method = CA.SmoothMinimumBlending(),
    ),
    turbconv_model = CA.PrognosticEDMFX(;
        n_updrafts = 1,
        prognostic_tke = true,
        area_fraction = FT(1.0e-3),
    ),
    # surface
    flux_scheme = CA.SurfaceConditions.MoninObukhov(;
        z0 = FT(surface_roughness()),
    ),
    temperature = SOCRATESSurfaceTemperature(),
    # numerics
    diff_mode = CA.Implicit(),
    hyperdiff = nothing,
    edmfx_sgsflux_upwinding = :first_order,
    kwargs...,
) where {FT <: AbstractFloat}
    return CA.AtmosModel(;
        external_forcing,
        insolation,
        microphysics_model,
        cloud_model,
        microphysics_tendency_timestepping,
        tracer_nonnegativity_method,
        sgs_quadrature,
        terminal_velocity_liquid,
        terminal_velocity_ice,
        terminal_velocity_rain,
        terminal_velocity_snow,
        radiation_mode,
        edmfx_model,
        turbconv_model,
        flux_scheme,
        temperature,
        diff_mode,
        hyperdiff,
        edmfx_sgsflux_upwinding,
        kwargs...,
    )
end

"""
    socrates_ode_config(FT; kwargs...)

The IMEX time-stepping algorithm for a SOCRATES case. Every argument of
`ClimaAtmos.ode_configuration` is exposed.
"""
socrates_ode_config(
    ::Type{FT};
    ode_algo::AbstractString = "ARS222",
    update_jacobian_every::AbstractString = "stage",
    max_newton_iters::Int = 1,
    use_krylov_method::Bool = false,
    use_dynamic_krylov_rtol::Bool = false,
    eisenstat_walker_forcing_alpha::Real = 2.0,
    krylov_rtol::Real = 0.1,
    use_newton_rtol::Bool = false,
    newton_rtol::Real = 1.0e-5,
    jvp_step_adjustment::Real = 1.0,
) where {FT <: AbstractFloat} = CA.ode_configuration(
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

# --- Simulation ------------------------------------------------------------- #

"""
    socrates_simulation(FT, case; output_dir, kwargs...)

Assemble an `AtmosSimulation` for `case`. Case values default to the shipped case
data and are all overridable; `kwargs...` are forwarded to `AtmosSimulation`.

`t_start` is fixed at 0: the SSCF forcing axis is relative to the case start, so
a nonzero start would silently offset the forcing from the state.
"""
function socrates_simulation(
    ::Type{FT},
    c::SocratesCase;
    output_dir::AbstractString,
    params = nothing,
    grid = socrates_grid(FT, c),
    t_end::Real = SOCRATES.t_end(c),
    dt::Real = 10,
    dt_rad = "10mins",
    forcing_source = c,
    scalar_nudge_timescale::Real = SOCRATES.scalar_nudge_timescale(),
    wind_nudge_timescale::Real = SOCRATES.wind_nudge_timescale(forcing_label(c)),
    diagnostics = socrates_diagnostics(; n_levels = length(socrates_z(grid))),
    insolation = SOCRATESInsolation(FT, c.flight_number, c.forcing_type),
    ode_config = socrates_ode_config(FT),
    jacobian = CA.ManualSparseJacobian(; approximate_solve_iters = 2),
    job_id::AbstractString = case_name(c),
    verbose::Bool = true,
    model_kwargs = (;),
    kwargs...,
) where {FT <: AbstractFloat}
    validate(c)
    z = socrates_z(grid)
    toml_dict = socrates_toml_dict(FT, c; params)
    atmos_params = socrates_params(toml_dict)
    thermo_params = CA.Parameters.thermodynamics_params(atmos_params)

    forcing = SOCRATESForcing(
        FT,
        forcing_source;
        scalar_nudge_timescale,
        wind_nudge_timescale,
    )
    setup = SOCRATESSetup(FT, c.flight_number, c.forcing_type, z; thermo_params)

    return CA.AtmosSimulation{FT}(;
        model = socrates_model(
            FT,
            atmos_params;
            external_forcing = forcing,
            insolation,
            terminal_velocity_modes(FT, toml_dict)...,
            model_kwargs...,
        ),
        params = atmos_params,
        grid,
        setup,
        dt,
        t_start = 0,
        t_end,
        start_date = les_start_datetime(c),
        ode_config,
        jacobian,
        diagnostics,
        callback_kwargs = (; dt_rad),
        job_id,
        output_dir,
        verbose,
        kwargs...,
    )
end

"""
    run_case(case; FT, output_dir, kwargs...) -> output_dir

Build and solve one SOCRATES case, returning the directory its diagnostics were
written to.
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
        "$(case_name(c)) did not run to completion: ClimaAtmos returned \
         `$(result.ret_code)`.",
    )
    return sim.output_dir
end

"""
    run_cases(cases; FT, output_dir, kwargs...) -> Vector

Run each case into its own subdirectory of `output_dir`. A case that fails stops
the sweep.
"""
run_cases(
    cases::AbstractVector{<:SocratesCase};
    FT::Type{<:AbstractFloat} = Float64,
    output_dir::AbstractString,
    kwargs...,
) = [
    run_case(c; FT, output_dir = joinpath(output_dir, case_name(c)), kwargs...)
    for c in cases
]
