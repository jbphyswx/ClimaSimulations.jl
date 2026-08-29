"""
    model.jl

The `AtmosModel` and `AtmosSimulation` for a CloudBench site.

The case identity, the sounding, the forcing, the initial condition, the insolation, the
surface condition and the ClimaParams overrides all come from `SwirlLMCloudBench`. What is
assembled here is the column around them.
"""

"""Duration [s] of a CloudBench LES, from `Simulation.CLOUDBENCH_LES_GRID`."""
cloudbench_t_end() = S.CLOUDBENCH_LES_GRID.duration_days * 86400.0

"""
Interval [s] between stored samples in `data.zarr`, so a run written on this cadence lines up
with the reference without resampling in time. `data.zarr` covers the last simulated day only.
"""
cloudbench_output_interval() = S.CLOUDBENCH_OUTPUT_INTERVAL

"""
Condensate, precipitation and the top-of-atmosphere radiative fluxes, all-sky and clear-sky.
A default for [`cloudbench_diagnostics`](@ref); pass `short_names` for any other set.
"""
const default_diagnostic_vars = (
    "clw", "husra", "cl", "lwp", "rwp", "pr",
    "rlut", "rsut", "rsdt", "rlutcs", "rsutcs",
    "hfls", "hfss",
)

"""The state a column is diagnosed from. A default for [`cloudbench_diagnostics`](@ref)."""
const STATE_VARS = ("ta", "hus", "ua", "va", "wa", "pfull", "rhoa")

"""
    cloudbench_diagnostics(short_names; period_seconds, reduction, n_levels, kwargs...)

A `DiagnosticsConfig` writing `short_names` on the model's own levels, on the reference
output cadence. `kwargs...` are forwarded to `DiagnosticsConfig`.

Names are checked against `ClimaAtmos`'s registry here so that a typo fails before the model
is built rather than part-way into a run.
"""
function cloudbench_diagnostics(
    short_names = (default_diagnostic_vars..., STATE_VARS...);
    period_seconds::Real = cloudbench_output_interval(),
    reduction::AbstractString = "average",
    n_levels::Int,
    kwargs...,
)
    isempty(short_names) && error("A run with no diagnostics writes no output.")
    unknown = [n for n in short_names if !haskey(CA.Diagnostics.ALL_DIAGNOSTICS, String(n))]
    isempty(unknown) ||
        error("Not ClimaAtmos diagnostics: $(join(unknown, ", ")).")
    return CA.DiagnosticsConfig(;
        default = false,
        additional = ((;
            short_name = collect(String.(short_names)),
            period = "$(Int(period_seconds))secs",
            reduction,
        ),),
        interpolation_num_points = (2, 2, n_levels),
        output_at_levels = true,
        kwargs...,
    )
end

"""
    cloudbench_sponge(FT; depth, domain_top, α_w)

A Rayleigh sponge on vertical velocity over the top `depth` metres, which is what the
reference damps (`Simulation.CLOUDBENCH_SPONGE_DEPTH`).

`zd` is measured from the ground, so it is the domain top less the sponge depth. Only `α_w`
is set: the reference damps `w` alone.

ClimaAtmos's profile is `α sin²(π(z-zd)/(2(zmax-zd)))`, the same shape swirl-lm applies, so
`α_w` is the reference's damping rate. swirl-lm's `β` is dimensionless — `1/a_coeff`, applied
as `β/dt` — so the rate is `1/(a_coeff dt)`, which is 0.25 s^-1 at its `a_coeff` default of 20
and the LES timestep of 0.2 s.
"""
function cloudbench_sponge(
    ::Type{FT};
    depth::Real = S.CLOUDBENCH_SPONGE_DEPTH,
    domain_top::Real = S.CLOUDBENCH_LES_GRID.lz,
    α_w::Real = 1 / (20 * S.CLOUDBENCH_LES_GRID.dt),
) where {FT <: AbstractFloat}
    0 < depth < domain_top ||
        error("A sponge depth of $depth m does not fit inside a $domain_top m domain.")
    return CA.AtmosSponge(;
        rayleigh_sponge = CA.RayleighSponge{FT}(;
            zd = FT(domain_top - depth),
            α_uₕ = FT(0),
            α_w = FT(α_w),
            α_tracer = FT(0),
        ),
    )
end

"""
    cloudbench_model(FT, params; external_forcing, insolation, kwargs...)

The `AtmosModel` for a CloudBench site: prognostic EDMFX with one updraft, 1-moment warm-rain
microphysics, interactive RRTMGP radiation, and a Rayleigh sponge over the top kilometre.

`kwargs...` are forwarded to `AtmosModel` and override any default below, so any component
can be replaced without editing this function.

The microphysics is the one place this cannot reproduce the reference exactly. CloudBench
runs single-moment warm rain with **instantaneous** condensation
(`Simulation.CLOUDBENCH_MICROPHYSICS`) and ClimaAtmos carries no equilibrium 1-moment model:
`NonEquilibriumMicrophysics1M` relaxes condensation over a finite timescale, and
`EquilibriumMicrophysics0M` removes precipitation instantly rather than carrying rain. The
1-moment model is the closer of the two and is the default; its condensation timescale is a
ClimaParams entry, so shortening it is a parameter change rather than a code change.
"""
function cloudbench_model(
    ::Type{FT},
    params;
    external_forcing,
    insolation,
    microphysics_model = CA.NonEquilibriumMicrophysics1M(;
        n_substeps = 3,
        n_substeps_quad = 2,
    ),
    cloud_model = CA.QuadratureCloud(),
    microphysics_tendency_timestepping = CA.Implicit(),
    sgs_quadrature = CA.SGSQuadrature(
        FT;
        quadrature_order = 3,
        distribution = CA.GaussianSGS(),
        T_min = FT(CA.Parameters.T_min_sgs(params)),
        q_max = FT(CA.Parameters.q_max_sgs(params)),
    ),
    radiation_mode = CA.RRTMGPInterface.AllSkyRadiationWithClearSkyDiagnostics(),
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
    sponge = cloudbench_sponge(FT),
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
        sgs_quadrature,
        radiation_mode,
        edmfx_model,
        turbconv_model,
        sponge,
        diff_mode,
        hyperdiff,
        edmfx_sgsflux_upwinding,
        kwargs...,
    )
end

"""
    cloudbench_ode_config(FT; kwargs...)

The IMEX time-stepping algorithm for a CloudBench case. Every argument of
`ClimaAtmos.ode_configuration` is exposed.
"""
cloudbench_ode_config(
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

"""A run directory name identifying the site, month and experiment."""
cloudbench_job_id(inst::S.CloudBenchInstance) =
    "cloudbench_$(inst.site_id)_$(inst.month)_$(Symbol(inst.experiment))"

cloudbench_job_id(sim::S.CloudBenchSimulation) =
    cloudbench_job_id(S.cloudbench_instance(sim))

"""
    cloudbench_start_date(case; year)

The first of the case's own month. `year` and the day label the run: CloudBench prescribes a
diurnally averaged TOA insolation and zenith angle per case, which the setup's insolation
model reads in place of the date.
"""
cloudbench_start_date(inst::S.CloudBenchInstance; year::Integer = 0000) =
    Dates.DateTime(year, inst.month, 1)

cloudbench_start_date(sim::S.CloudBenchSimulation; kwargs...) =
    cloudbench_start_date(S.cloudbench_instance(sim); kwargs...)

"""
    cloudbench_simulation(FT, case; output_dir, kwargs...)

Assemble an `AtmosSimulation` for a CloudBench `case`, given as a
`Simulation.CloudBenchInstance` or `Simulation.CloudBenchSimulation`.

The setup, the parameters and the radiation cadence come from `SwirlLMCloudBench`, so this
case's experiment CO₂, its SST, its zenith angle and its TOA insolation are the case's own.
`kwargs...` are forwarded to `AtmosSimulation`.

`t_start` is 0: the sounding is steady and its forcing carries no time axis, so a nonzero
start would mean nothing.
"""
function cloudbench_simulation(
    ::Type{FT},
    case;
    output_dir::AbstractString,
    grid = cloudbench_grid(FT),
    params = SW.ClimaAtmos_SwirlLMCloudBench_params(case, FT),
    setup = SW.ClimaAtmosSwirlLMCloudBenchSetup(case; FT),
    t_end::Real = cloudbench_t_end(),
    dt::Real = 20,
    diagnostics = cloudbench_diagnostics(; n_levels = length(cloudbench_z(grid))),
    ode_config = cloudbench_ode_config(FT),
    jacobian = CA.ManualSparseJacobian(; approximate_solve_iters = 2),
    start_date = cloudbench_start_date(case),
    job_id::AbstractString = cloudbench_job_id(case),
    callback_kwargs = SW.ClimaAtmos_SwirlLMCloudBench_callback_kwargs(),
    verbose::Bool = true,
    model_kwargs = (;),
    kwargs...,
) where {FT <: AbstractFloat}
    return CA.AtmosSimulation{FT}(;
        model = cloudbench_model(
            FT,
            params;
            external_forcing = CA.Setups.external_forcing(setup, FT),
            insolation = CA.Setups.insolation_model(setup),
            model_kwargs...,
        ),
        params,
        grid,
        setup,
        dt,
        t_start = 0,
        t_end,
        start_date,
        ode_config,
        jacobian,
        diagnostics,
        callback_kwargs,
        job_id,
        output_dir,
        verbose,
        kwargs...,
    )
end

"""
    run_case(case; FT, output_dir, kwargs...) -> output_dir

Build and solve one CloudBench case, returning the directory its diagnostics were written to.
"""
function run_case(
    case;
    FT::Type{<:AbstractFloat} = Float64,
    output_dir::AbstractString,
    kwargs...,
)
    sim = cloudbench_simulation(FT, case; output_dir, kwargs...)
    result = CA.solve_atmos!(sim)
    result.ret_code === :success || error(
        "$(cloudbench_job_id(case)) did not run to completion: ClimaAtmos returned \
         `$(result.ret_code)`.",
    )
    return output_dir
end