"""
    model.jl

The one model component this case owns — its frozen insolation — and the
assembly of an `AtmosSimulation` for one AYiL day.

Every default here is overridable: `mosaic_model` and `mosaic_simulation` forward
`kwargs...` to `AtmosModel` and `AtmosSimulation`.
"""

# --- Insolation ------------------------------------------------------------- #

"""
    MOSAiCInsolation{FT}(cos_zenith, toa_flux)

Insolation held at a fixed zenith angle, as the reference runs did
(`lcnstzenithtime = .true.`, `cnstzenithtime = 11` on all 190 days).

`toa_flux` is the flux perpendicular to the beam, so RRTMGP's downwelling
shortwave is `toa_flux * cos_zenith`.
"""
struct MOSAiCInsolation{FT} <: CA.AbstractInsolation
    cos_zenith::FT
    toa_flux::FT
end

"""
    MOSAiCInsolation(FT, case; reference_time, latitude, longitude, insolation_params)

The orbital insolation at the case's reference time and drift position.

Where the sun is below the horizon the pair is `(eps(FT), 0)`: no incoming flux,
with the positive zenith cosine RRTMGP requires. At 77-90 °N a large share of the
190 days are polar night, so that is an ordinary case here, not an edge one.
"""
function MOSAiCInsolation(
    ::Type{FT},
    c::MOSAiCAYiLCase;
    root = data_root(),
    forcing = read_scm_in(c.date; root),
    reference_time::Dates.DateTime = reference_datetime(c),
    latitude::Real = forcing.surface.trajectory_latitude,
    longitude::Real = forcing.surface.trajectory_longitude,
    insolation_params = Insolation.Parameters.InsolationParameters(FT),
) where {FT <: AbstractFloat}
    (; S, μ) =
        Insolation.insolation(reference_time, latitude, longitude, insolation_params)
    return μ > 0 ? MOSAiCInsolation{FT}(FT(μ), FT(S)) :
           MOSAiCInsolation{FT}(eps(FT), zero(FT))
end

function CA.set_insolation_variables!(Y, p, t, insolation::MOSAiCInsolation)
    (; rrtmgp_solver) = p.radiation
    CA.RRTMGP.cos_zenith(rrtmgp_solver) .= insolation.cos_zenith
    CA.RRTMGP.toa_flux(rrtmgp_solver) .= insolation.toa_flux
    return nothing
end

"""
    reference_datetime(case; hour = 11)

The case's reference time: `cnstzenithtime` UTC on its date, the hour the frozen
zenith angle and the initializing radiosonde both refer to.
"""
reference_datetime(c::MOSAiCAYiLCase; hour::Integer = 11) =
    Dates.DateTime(Dates.Date(c)) + Dates.Hour(hour)

# --- Diagnostics ------------------------------------------------------------ #

"""
The variables written by default: the condensate profiles and their column
integrals, the thermodynamic profiles, and the radiative fluxes, which are the
quantities the DALES output carries a counterpart for (`docs/design.md` §9, §11).
"""
const DEFAULT_DIAGNOSTIC_VARS = (
    "clw", "cli", "husra", "hussn",
    "ql_all", "qi_all",
    "lwp", "iwp", "rwp", "swp",
    "ta", "hus",
    # density, because a water path is integrated from `ρq` on the comparison grid
    # rather than differenced against a stored path (`docs/design.md` §11)
    "rhoa",
    "rlut", "rsut", "rlds", "rsds", "rlus", "rsus",
    # face flux profiles, all-sky and clear-sky: the reference writes exactly these
    # (`lwd`/`lwu`/`swd`/`swu` and the clear-air `*ca` pair), so radiation is
    # compared profile against profile with no heating rate reconstructed
    "rld", "rlu", "rsd", "rsu",
    "rldcs", "rlucs", "rsdcs", "rsucs",
    # the surface scheme's own output, against `wthls`/`wqts` (`surface_heat_fluxes`)
    "hfss", "hfls",
)

"""
    mosaic_diagnostics(short_names; period_seconds, reduction, n_levels, kwargs...)

A `DiagnosticsConfig` writing `short_names` on the model's own levels.

`period_seconds` defaults to the reference's own output interval, which writes
300 s averages sampled every 60 s (`docs/design.md` §11).
"""
mosaic_diagnostics(
    short_names = DEFAULT_DIAGNOSTIC_VARS;
    period_seconds::Real = 300,
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
    mosaic_model(FT, params; external_forcing, insolation, scm_coriolis, temperature, flux_scheme, surface_albedo, kwargs...)

The `AtmosModel` for an AYiL day: prognostic EDMFX with one updraft,
non-equilibrium 1-moment microphysics with a quadrature cloud, interactive
RRTMGP radiation, and Monin-Obukhov surface fluxes over a prescribed skin
temperature.

The reference ran Seifert-Beheng two-moment microphysics, a different scheme;
`docs/design.md` §13 records that and the other differences.

Every keyword is a whole model component, so each concept has exactly one way to
be set: to change the number of updrafts, pass a `turbconv_model`, not a separate
`n_updrafts` that would silently contend with it. The defaults are default
*arguments*, so a component that is passed in is never also constructed. Any
`AtmosModel` field not named here can be set through `kwargs...`.

The six components with no default are the ones that need case data: the caller
has already read it, and guessing it here would mean reading the file again.
"""
function mosaic_model(
    ::Type{FT},
    params;
    # from the case's data
    external_forcing,
    insolation,
    scm_coriolis,
    temperature,
    flux_scheme,
    surface_albedo,
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
    # numerics
    diff_mode = CA.Implicit(),
    hyperdiff = nothing,
    edmfx_sgsflux_upwinding = :first_order,
    kwargs...,
) where {FT <: AbstractFloat}
    return CA.AtmosModel(;
        external_forcing,
        insolation,
        scm_coriolis,
        temperature,
        flux_scheme,
        surface_albedo,
        microphysics_model,
        cloud_model,
        microphysics_tendency_timestepping,
        tracer_nonnegativity_method,
        sgs_quadrature,
        radiation_mode,
        edmfx_model,
        turbconv_model,
        diff_mode,
        hyperdiff,
        edmfx_sgsflux_upwinding,
        kwargs...,
    )
end

# --- Simulation ------------------------------------------------------------- #

"""
    mosaic_simulation(FT, case; output_dir, kwargs...)

Assemble an `AtmosSimulation` for one AYiL day.

`scm_in` is read once here and shared by the forcing, the initial condition, the
insolation and the Coriolis term, so a case costs one open of the file rather
than one per consumer.

`grid` is the reference's own full column by default. To avoid problems from anoamlous ice
(see documentation) manually consruct a grid that stops at [`best_simulation_top`](@ref), e.g.:

```julia
faces = truncate_faces_to_top(native_faces(c), best_simulation_top(c))
sim = mosaic_simulation(FT, c; output_dir, grid = mosaic_grid(FT, c; faces))
```

`t_start` is fixed at 0: the forcing and `t_end` are both relative to the run
start.
"""
function mosaic_simulation(
    ::Type{FT},
    c::MOSAiCAYiLCase;
    output_dir::AbstractString,
    root = data_root(),
    time_index::Int = 1,
    forcing_data = read_scm_in(c.date; root, time_index),
    params = mosaic_params(FT, c),
    grid = mosaic_grid(FT, c; root),
    t_end::Real = MOSAiC_AYiL.t_end(c; root),
    dt::Real = 10,
    dt_rad = "10mins",
    z0::Real = forcing_data.surface.z0_momentum,
    albedo::Real = forcing_data.surface.albedo,
    ode_config = CA.ode_configuration(
        FT, "ARS222", "stage", 1, false, false, 2.0, 0.1, false, 1.0e-5, 1.0,
    ),
    jacobian = CA.ManualSparseJacobian(; approximate_solve_iters = 2),
    diagnostics = mosaic_diagnostics(; n_levels = length(mosaic_z(grid))),
    job_id::AbstractString = case_name(c),
    verbose::Bool = true,
    model_kwargs = (;),
    kwargs...,
) where {FT <: AbstractFloat}
    z = mosaic_z(grid)
    return CA.AtmosSimulation{FT}(;
        model = mosaic_model(
            FT,
            params;
            external_forcing = MOSAiCForcing(
                FT, c; root, time_index, forcing = forcing_data,
            ),
            insolation = MOSAiCInsolation(FT, c; root, forcing = forcing_data),
            scm_coriolis = mosaic_scm_coriolis(
                FT, c; params, root, time_index, forcing = forcing_data,
            ),
            temperature = CA.SurfaceConditions.AnalyticTemperature(
                Returns(FT(surface_temperature(forcing_data))),
            ),
            flux_scheme = CA.SurfaceConditions.MoninObukhov(; z0 = FT(z0)),
            surface_albedo = CA.ConstantAlbedo{FT}(; α = FT(albedo)),
            model_kwargs...,
        ),
        params,
        grid,
        setup = MOSAiCSetup(FT, c; root, time_index, forcing = forcing_data),
        dt,
        t_start = 0,
        t_end,
        start_date = reference_datetime(c),
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

Build and solve one AYiL day, returning the directory its diagnostics were
written to.
"""
function run_case(
    c::MOSAiCAYiLCase;
    FT::Type{<:AbstractFloat} = Float64,
    output_dir::AbstractString,
    kwargs...,
)
    sim = mosaic_simulation(FT, c; output_dir, kwargs...)
    result = CA.solve_atmos!(sim)
    result.ret_code === :success || error(
        "$(case_name(c)) did not run to completion: ClimaAtmos returned \
         `$(result.ret_code)`.",
    )
    return sim.output_dir
end
