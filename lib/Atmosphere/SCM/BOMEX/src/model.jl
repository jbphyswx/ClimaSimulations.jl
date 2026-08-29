"""
Condensate, precipitation and the surface fluxes. A default for [`bomex_diagnostics`](@ref);
pass `short_names` for any other set.

BOMEX carries no radiation, so no radiative flux is written — asking for one would ask the model
for a quantity it does not compute.
"""
const default_diagnostic_vars =
    ("clw", "cli", "husra", "hussn", "cl", "lwp", "rwp", "pr", "hfls", "hfss")

"""The state a column is diagnosed from. A default for [`bomex_diagnostics`](@ref)."""
const STATE_VARS = ("ta", "thetaa", "hus", "ua", "va", "wa", "pfull", "rhoa")

"""
    bomex_diagnostics(short_names; period_seconds, reduction, n_levels, kwargs...)

A `DiagnosticsConfig` writing `short_names` on the model's own levels. `kwargs...` are
forwarded to `DiagnosticsConfig`.

Names are checked against `ClimaAtmos`'s registry here, so a typo fails before the model is
built rather than part-way into a run.
"""
function bomex_diagnostics(
    short_names = (default_diagnostic_vars..., STATE_VARS...);
    period_seconds::Real = 600,
    reduction::AbstractString = "average",
    n_levels::Int,
    kwargs...,
)
    isempty(short_names) && error("A run with no diagnostics writes no output.")
    unknown =
        [n for n in short_names if !haskey(CA.Diagnostics.ALL_DIAGNOSTICS, String(n))]
    isempty(unknown) || error("Not ClimaAtmos diagnostics: $(join(unknown, ", ")).")
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
    bomex_model(FT, params, c; kwargs...)

The `AtmosModel` for `c`: prognostic EDMFX with single-moment microphysics, driven by the
case's prescribed subsidence, large-scale advection and Coriolis forcing.

No radiation: BOMEX prescribes its radiative cooling as part of the large-scale temperature
tendency, so an interactive scheme would double-count it.
"""
function bomex_model(
    ::Type{FT},
    params,
    c::BOMEXCase;
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
        prognostic_tke = c.prognostic_tke,
        area_fraction = FT(1.0e-3),
    ),
    diff_mode = CA.Implicit(),
    hyperdiff = nothing,
    edmfx_sgsflux_upwinding = :first_order,
    kwargs...,
) where {FT <: AbstractFloat}
    return CA.AtmosModel(;
        subsidence = CA.get_subsidence_model(FT; setup_type = c),
        ls_adv = CA.get_large_scale_advection_model(FT; setup_type = c),
        scm_coriolis = CA.get_scm_coriolis(FT; setup_type = c),
        microphysics_model,
        cloud_model,
        microphysics_tendency_timestepping,
        sgs_quadrature,
        edmfx_model,
        turbconv_model,
        diff_mode,
        hyperdiff,
        edmfx_sgsflux_upwinding,
        kwargs...,
    )
end

"""
    bomex_ode_config(FT; kwargs...)

The IMEX time-stepping algorithm for BOMEX. Every argument of
`ClimaAtmos.ode_configuration` is exposed.
"""
bomex_ode_config(
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

"""A run directory name for the case."""
bomex_job_id(c::BOMEXCase) = case_name(c)

"""
    bomex_simulation(FT, c; output_dir, kwargs...)

The `AtmosSimulation` for `c`.

`t_start` is 0: the case's forcing is steady and carries no time axis, so a nonzero start would
mean nothing.
"""
function bomex_simulation(
    ::Type{FT},
    c::BOMEXCase;
    output_dir::AbstractString,
    params = bomex_params(FT),
    grid = bomex_grid(FT, c),
    t_end::Real = BOMEX.t_end(c),
    dt::Real = 120,
    diagnostics = bomex_diagnostics(; n_levels = length(bomex_z(grid))),
    ode_config = bomex_ode_config(FT),
    jacobian = CA.ManualSparseJacobian(; approximate_solve_iters = 2),
    job_id::AbstractString = bomex_job_id(c),
    verbose::Bool = true,
    model_kwargs = (;),
    kwargs...,
) where {FT <: AbstractFloat}
    return CA.AtmosSimulation{FT}(;
        model = bomex_model(FT, params, c; model_kwargs...),
        params,
        grid,
        setup = c,
        dt,
        t_start = 0,
        t_end,
        ode_config,
        jacobian,
        diagnostics,
        job_id,
        output_dir,
        verbose,
        kwargs...,
    )
end

"""
    run_case(c; FT, output_dir, kwargs...) -> output_dir

Build and solve BOMEX, returning the directory its diagnostics were written to.
"""
function run_case(
    c::BOMEXCase;
    FT::Type{<:AbstractFloat} = Float64,
    output_dir::AbstractString,
    kwargs...,
)
    sim = bomex_simulation(FT, c; output_dir, kwargs...)
    result = CA.solve_atmos!(sim)
    result.ret_code === :success || error(
        "$(bomex_job_id(c)) did not run to completion: ClimaAtmos returned \
         `$(result.ret_code)`.",
    )
    return sim.output_dir
end