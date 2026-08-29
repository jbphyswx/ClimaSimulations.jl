"""
    calibration.jl

The EKI interface for BOMEX: sampled parameters into a forward run, run output back out as a `G`
ensemble, and the prior and observations the inversion needs.

Nothing here is reachable from `BOMEX`. A member's parameters are built by layering the sampled
file over the case's own `configs/BOMEX.toml` through `BOMEX.bomex_params`, so the forward model
needs no knowledge that a calibration exists.
"""

# --- Sampled parameters into a run ---------------------------------------------------------- #

"""
    member_params(parameter_file; FT)

`ClimaAtmosParameters` for one ensemble member: the case's own overrides with the sampled values
layered on top.
"""
function member_params(
    parameter_file::AbstractString;
    FT::Type{<:AbstractFloat} = Float64,
)
    isfile(parameter_file) || error("No sampled parameter file at $parameter_file.")
    return BOMEX.bomex_params(FT; params = parameter_file)
end

# --- The prior ------------------------------------------------------------------------------ #

"""
Parameter => `(mean, lower, upper, unconstrained_σ)`.

A starting point, not a prescription: pass your own `spec` to [`default_prior`](@ref).
"""
const DEFAULT_PRIOR_SPEC = (
    condensation_evaporation_timescale = (1.0e2, 1.0e0, 1.0e8, 3.0),
    entr_coeff = (1.0e-1, 1.0e-3, 1.0e1, 2.0),
    detr_massflux_vertdiv_coeff = (3.0e-1, 1.0e-3, 1.0e1, 2.0),
)

"""
    default_prior(spec = DEFAULT_PRIOR_SPEC)

A `ParameterDistribution` over `spec`, each entry bounded and normal in the unconstrained space.
"""
function default_prior(spec = DEFAULT_PRIOR_SPEC)
    isempty(spec) && error("The prior needs at least one parameter.")
    distributions = map(collect(pairs(spec))) do (name, (mean, lower, upper, σ))
        lower <= mean <= upper ||
            error("Prior mean $mean for `$name` is outside bounds ($lower, $upper).")
        σ > 0 || error("Prior `unconstrained_σ` for `$name` must be positive, got $σ.")
        constraint = EKP.ParameterDistributions.bounded(lower, upper)
        return EKP.ParameterDistributions.ParameterDistribution(
            EKP.ParameterDistributions.Parameterized(
                EKP.ParameterDistributions.Normal(
                    constraint.constrained_to_unconstrained(mean), σ,
                ),
            ),
            constraint,
            String(name),
        )
    end
    return EKP.ParameterDistributions.combine_distributions(distributions)
end

"""The parameter names [`default_prior`](@ref) samples, in order."""
prior_names(spec = DEFAULT_PRIOR_SPEC) = collect(String.(keys(spec)))

# --- The interface -------------------------------------------------------------------------- #

"""
    BOMEXInterface(; reference, output_dir, case, vars, float_type, grid, run_kwargs)

`ClimaCalibrate` model interface for a BOMEX EKI calibration.

`reference` is required and has no default: this case ships none, so what a member is scored
against is the caller's to supply. See [`read_reference`](@ref).

`run_kwargs` are forwarded to `BOMEX.run_case`, so `t_end`, `dt` and the model components are
the caller's to set.
"""
struct BOMEXInterface{FT <: AbstractFloat, C, R, G, K} <:
       ClimaCalibrate.AbstractModelInterface
    case::C
    reference::R
    output_dir::String
    vars::Vector{String}
    float_type::Type{FT}
    grid::G
    run_kwargs::K
    startup_waves::Int
    startup_wave_pause::Float64
end

function BOMEXInterface(;
    reference,
    output_dir::AbstractString,
    float_type::Type{<:AbstractFloat} = Float64,
    case = BOMEX.case(float_type),
    vars = collect(String, default_calibration_vars),
    grid = BOMEX.bomex_grid(float_type, case),
    run_kwargs = (;),
    startup_waves::Integer = 0,
    startup_wave_pause::Real = 0,
)
    isempty(vars) && error("BOMEXInterface needs at least one scored variable.")
    haskey(run_kwargs, :params) &&
        error("`params` is set per member from the sampled file; do not pass it.")
    any(k -> haskey(run_kwargs, k), (:grid, :faces)) &&
        error("Pass the vertical grid as `grid`, not `run_kwargs`.")
    output_dir = abspath(output_dir)
    mkpath(output_dir)
    return BOMEXInterface(
        case,
        reference,
        output_dir,
        collect(String, vars),
        float_type,
        grid,
        run_kwargs,
        Int(startup_waves),
        Float64(startup_wave_pause),
    )
end

"""The heights a member is scored on: the model grid's own centres inside the compared region."""
comparison_levels(interface::BOMEXInterface) = scored_levels(
    BOMEX.bomex_z(interface.grid),
    default_z_bounds(interface.case),
)

member_output_dir(interface::BOMEXInterface, iteration, member) = joinpath(
    ClimaCalibrate.path_to_ensemble_member(interface.output_dir, iteration, member),
    BOMEX.case_name(interface.case),
)

# --- Observations --------------------------------------------------------------------------- #

"""
    observation(interface)

The `EKP.Observation` the inversion is against: each scored variable on the comparison levels,
stacked in `interface.vars` order.
"""
function observation(
    interface::BOMEXInterface,
    ::Type{FT} = Float64;
    uncertainty_fraction::Real = 0.1,
    floor_fraction::Real = 1.0e-3,
) where {FT <: AbstractFloat}
    uncertainty_fraction > 0 ||
        error("`uncertainty_fraction` must be positive, got $uncertainty_fraction.")
    levels = comparison_levels(interface)
    truth = reference_values(interface.reference, interface.vars, levels)
    y = FT[]
    diagonal = FT[]
    for name in interface.vars
        values = truth[name]
        scale = maximum(abs, values)
        scale > 0 || error(
            "Reference `$name` is identically zero over the compared levels, so it carries no \
             weight; drop it from `vars`.",
        )
        for v in values
            push!(y, FT(v))
            # a floor, so a level where the reference is ~0 does not get an infinite weight
            σ = uncertainty_fraction * max(abs(v), floor_fraction * scale)
            push!(diagonal, FT(σ^2))
        end
    end
    return EKP.Observation(
        Dict(
            "samples" => y,
            "covariances" => LinearAlgebra.Diagonal(diagonal),
            "names" => BOMEX.case_name(interface.case),
        ),
    )
end

"""
    build_ekp(interface, prior; ensemble_size, kwargs...)

The `EnsembleKalmanProcess` for this interface.
"""
function build_ekp(
    interface::BOMEXInterface,
    prior = default_prior();
    ensemble_size::Int,
    rng = Random.MersenneTwister(1234),
    verbose::Bool = true,
    kwargs...,
)
    ensemble_size >= 2 ||
        error("ensemble_size must be at least 2, got $ensemble_size.")
    return EKP.EnsembleKalmanProcess(
        EKP.construct_initial_ensemble(rng, prior, ensemble_size),
        observation(interface, interface.float_type),
        EKP.Inversion();
        rng,
        verbose,
        kwargs...,
    )
end

# --- The ClimaCalibrate methods -------------------------------------------------------------- #

"""
    run_member(interface, iteration, member)

Run one member, returning the directory its diagnostics were written to. The unit of work
[`ClimaCalibrate.Calibration.run_iteration`](@ref) schedules, so it has to be callable on a
worker.
"""
function run_member(interface::BOMEXInterface, iteration::Integer, member::Integer)
    parameter_file =
        ClimaCalibrate.parameter_path(interface.output_dir, iteration, member)
    return BOMEX.run_case(
        interface.case;
        FT = interface.float_type,
        output_dir = member_output_dir(interface, iteration, member),
        grid = interface.grid,
        params = member_params(parameter_file; FT = interface.float_type),
        verbose = false,
        interface.run_kwargs...,
    )
end

ClimaCalibrate.forward_model(interface::BOMEXInterface, iteration, member) =
    (run_member(interface, iteration, member); nothing)

"""
    member_g_column(interface, iteration, member)

One member's contribution to `G`: the same reduction [`observation`](@ref) applies, taken over
that member's run output, stacked in `interface.vars` order.
"""
function member_g_column(
    interface::BOMEXInterface,
    iteration,
    member,
    ::Type{FT} = Float64,
) where {FT <: AbstractFloat}
    levels = comparison_levels(interface)
    values = model_values(
        member_output_dir(interface, iteration, member),
        interface.vars,
        levels;
        window = score_window(interface.case),
    )
    return FT[v for name in interface.vars for v in values[name]]
end

function ClimaCalibrate.observation_map(interface::BOMEXInterface, iteration)
    ekp = ClimaCalibrate.load_ekp_struct(interface.output_dir)
    n_members = EKP.get_N_ens(ekp)
    g_ensemble = fill(NaN, length(EKP.get_obs(ekp)), n_members)
    for member in 1:n_members
        try
            g_ensemble[:, member] = member_g_column(interface, iteration, member)
        catch e
            @warn "Member $member failed the observation map; leaving its column NaN" exception =
                (e, catch_backtrace())
        end
    end
    all(isnan, g_ensemble) &&
        error("Every member of iteration $iteration produced a NaN column.")
    return g_ensemble
end

ClimaCalibrate.model_interface_filepath(::BOMEXInterface) =
    abspath(joinpath(@__DIR__, "BOMEXCalibration.jl"))

ClimaCalibrate.experiment_dir(::BOMEXInterface) = abspath(joinpath(@__DIR__, "..", ".."))

# --- Flat parallelism ------------------------------------------------------------------------ #

member_marker_path(interface::BOMEXInterface, iteration, member) =
    joinpath(member_output_dir(interface, iteration, member), "member_completed")

member_completed(interface::BOMEXInterface, iteration, member) =
    isfile(member_marker_path(interface, iteration, member))

function mark_member_completed(interface::BOMEXInterface, iteration, member)
    path = member_marker_path(interface, iteration, member)
    mkpath(dirname(path))
    write(path, "completed")
    return nothing
end

"""
    ClimaCalibrate.Calibration.run_iteration(backend::WorkerBackend, interface::BOMEXInterface, ...)

Run one iteration over the worker pool.

BOMEX has a single case, so the flat `(member, case)` pool the multi-case siblings need
degenerates to one task per member — but it is still routed through `BOMEX.run_tasks` so the
startup-wave staggering and the resume-from-markers behaviour are the same as theirs.
"""
function ClimaCalibrate.Calibration.run_iteration(
    backend::ClimaCalibrate.WorkerBackend,
    interface::BOMEXInterface,
    iteration,
    ensemble_size,
    output_dir,
)
    tasks = [m for m in 1:ensemble_size if !member_completed(interface, iteration, m)]
    @info "Iteration $iteration: $(length(tasks))/$ensemble_size members to run" n_workers =
        length(backend.worker_pool.workers)

    if !isempty(tasks)
        executor = BOMEX.WorkerPoolExecutor(
            backend.worker_pool;
            empty_pool_timeout = backend.empty_pool_timeout,
            interface.startup_waves,
            interface.startup_wave_pause,
        )
        # one grid for every member, so every task carries the same affinity key
        keys = fill(BOMEX.topology_key(interface.grid), length(tasks))
        run_one = member -> run_member(interface, iteration, member)
        results = BOMEX.run_tasks(run_one, tasks, keys, executor)
        for (member, result) in zip(tasks, results)
            isnothing(result) || mark_member_completed(interface, iteration, member)
        end
    end

    failed = [m for m in 1:ensemble_size if !member_completed(interface, iteration, m)]
    for member in 1:ensemble_size
        member in failed ||
            ClimaCalibrate.write_model_completed(output_dir, iteration, member)
    end
    rate = length(failed) / ensemble_size
    rate > backend.failure_rate && error(
        "Iteration $iteration had a $(round(rate * 100; digits = 2))% member failure rate \
         (members $failed), exceeding the $(backend.failure_rate * 100)% threshold.",
    )
    isempty(failed) || @warn "Failed ensemble members: $failed"
    return nothing
end
