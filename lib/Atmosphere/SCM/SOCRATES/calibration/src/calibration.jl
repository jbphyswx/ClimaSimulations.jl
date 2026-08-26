"""
    calibration.jl

Ensemble Kalman Inversion (EKI) calibration interface, observation mapping, and driver for SOCRATES.
"""

using ClimaAnalysis: ClimaAnalysis
using ClimaCalibrate: ClimaCalibrate
using Dates: Dates
using Distributed: Distributed
using EnsembleKalmanProcesses: EnsembleKalmanProcesses as EKP
using JLD2: JLD2
using LinearAlgebra: LinearAlgebra
using Logging: Logging
using NaNStatistics: NaNStatistics
using Random: Random
using SOCRATES: SOCRATES
using Statistics: Statistics

# --- Observations Construction ------------------------------------------------------------- #

"""
    normalized_reference_series(case; transform, vars, window, bounds, float_type, grid, z_grid)

Stack reference observations row-wise, normalized by `pool_var`.
"""
function normalized_reference_series(
    c::SOCRATES.SOCRATESCase;
    transform::ScoreTransform = ScoreTransform(),
    vars = SOCRATES.SCORED_VARS,
    window = score_window(c),
    bounds = default_z_bounds(c),
    float_type::Type{<:AbstractFloat} = Float64,
    grid = SOCRATES.socrates_grid(float_type, c),
    z_grid = collect(Float64, SOCRATES.socrates_z(grid)),
)
    reference = les_outputvars(c; vars)
    rows = Vector{Matrix{Float64}}()
    pool_vars = Dict{String, Float64}()
    ranges = Vector{Pair{String, UnitRange{Int}}}()
    windowed = Dict{String, ClimaAnalysis.OutputVar}()
    offset = 0
    for name in vars
        var = _reference_on_grid(reference[name], z_grid, window, bounds)
        windowed[name] = var
        series = _series_matrix(var)
        mean_profile = vec(Statistics.mean(series; dims = 2))
        pv = pool_var(transform, name, mean_profile)
        pool_vars[name] = pv
        push!(rows, series ./ sqrt(pv))
        push!(ranges, name => (offset + 1):(offset + size(series, 1)))
        offset += size(series, 1)
    end
    return (; series = reduce(vcat, rows), pool_vars, ranges, windowed, vars)
end

function _reference_on_grid(var, z_grid, window, bounds)
    out = var
    if ClimaAnalysis.has_altitude(out)
        out = reference_on_levels(out, scored_levels(z_grid, bounds))
    end
    return ClimaAnalysis.window(
        out,
        ClimaAnalysis.time_name(out);
        left = first(window),
        right = last(window),
    )
end

_series_matrix(var) =
    ClimaAnalysis.has_altitude(var) ? Array{Float64}(var.data) :
    reshape(Array{Float64}(vec(var.data)), 1, :)

"""
    case_observation(case; transform, vars, rank, kwargs...)

Build `EKP.Observation` for one case with full cross-variable SVD + diagonal noise covariance.
"""
function case_observation(
    c::SOCRATES.SOCRATESCase;
    transform::ScoreTransform = ScoreTransform(),
    vars = SOCRATES.SCORED_VARS,
    rank = nothing,
    kwargs...,
)
    ref = normalized_reference_series(c; transform, vars, kwargs...)
    series = ref.series
    y = _nanmean_rows(series)
    diagonal = zeros(Float64, size(series, 1))
    for (name, range) in ref.ranges
        diagonal[range] .= uncertainty_diagonal(transform, name, view(series, range, :))
    end
    bad = findall(<=(0.0), diagonal)
    isempty(bad) || error(
        "Observation variance for $(SOCRATES.case_name(c)) is not positive at $(length(bad)) entries.",
    )
    Γ = EKP.SVDplusD(
        _low_rank_time_covariance(series, rank),
        LinearAlgebra.Diagonal(diagonal),
    )
    metadata = [
        _flat_metadata(ref.windowed[name], ref.pool_vars[name], scored_name(c, name))
        for name in vars
    ]
    return EKP.Observation(
        Dict(
            "samples" => y,
            "covariances" => Γ,
            "names" => SOCRATES.case_name(c),
            "metadata" => metadata,
        ),
    )
end

scored_name(c::SOCRATES.SOCRATESCase, var::AbstractString) =
    string(SOCRATES.case_name(c), "_", var)

function _flat_metadata(var, pool_var, short_name)
    tagged = ClimaAnalysis.OutputVar(
        merge(
            var.attributes,
            Dict{String, Any}(
                "pool_var" => pool_var,
                "units" => "1",
                "short_name" => short_name,
            ),
        ),
        var.dims,
        var.dim_attributes,
        var.data,
    )
    mean_var = ClimaAnalysis.average_time(tagged)
    return ClimaAnalysis.flatten(mean_var).metadata
end

_nanmean_rows(m) = [SOCRATES.nanmean(view(m, i, :)) for i in axes(m, 1)]

function _low_rank_time_covariance(series, rank)
    samples = replace(series, NaN => 0.0)
    return isnothing(rank) ? EKP.tsvd_cov_from_samples(samples; quiet = true) :
           EKP.tsvd_cov_from_samples(samples, rank; quiet = true)
end

"""
    observation_vector(cases; float_type, grids, kwargs...)

One `EKP.Observation` per case.
"""
function observation_vector(
    cases;
    float_type::Type{<:AbstractFloat} = Float64,
    grids = [SOCRATES.socrates_grid(float_type, c) for c in cases],
    kwargs...,
)
    length(grids) == length(cases) || error(
        "observation_vector needs one grid per case: got $(length(grids)) grids for $(length(cases)) cases.",
    )
    return [
        case_observation(c; float_type, grid, kwargs...) for
        (c, grid) in zip(cases, grids)
    ]
end

# --- Observation Map & GEnsembleBuilder ----------------------------------------------------- #

function model_scored_var(
    output_dir::AbstractString,
    name::AbstractString,
    c::SOCRATES.SOCRATESCase,
    pool_var::Real;
    window = score_window(c),
    bounds,
    reduction::AbstractString = "average",
    period::AbstractString = "10m",
)
    var = run_outputvars(output_dir, (name,); reduction, period)[name]
    levels = model_levels(var, bounds)
    mean_var = windowed_time_mean(restrict_to_levels(var, levels), window)
    data = similar(mean_var.data, Float64)
    data .= mean_var.data ./ sqrt(pool_var)
    dims = ClimaAnalysis.Var.OrderedDict{String, Vector{Float64}}(
        dim => collect(Float64, values) for (dim, values) in mean_var.dims
    )
    return ClimaAnalysis.OutputVar(
        merge(
            mean_var.attributes,
            Dict{String, Any}(
                "units" => "1",
                "short_name" => scored_name(c, name),
            ),
        ),
        dims,
        mean_var.dim_attributes,
        data,
    )
end

function pool_var_from_metadata(metadata)
    out = Dict{String, Float64}()
    for m in metadata
        name = get(m.attributes, "short_name", nothing)
        pv = get(m.attributes, "pool_var", nothing)
        isnothing(name) && error("An observation metadata entry has no `short_name`.")
        isnothing(pv) && error("Observation metadata for `$name` carries no `pool_var`.")
        out[name] = Float64(pv)
    end
    return out
end

function build_g_ensemble(interface, iteration::Integer)
    ekp = ClimaCalibrate.load_ekp_struct(interface.output_dir, iteration)
    builder = ClimaCalibrate.EnsembleBuilder.GEnsembleBuilder(ekp)
    metadata = ClimaCalibrate.get_metadata_for_nth_iteration(
        EKP.get_observation_series(ekp),
        EKP.get_N_iterations(ekp) + 1,
    )
    pool_vars = pool_var_from_metadata(metadata)

    for member in 1:EKP.get_N_ens(ekp)
        try
            for c in interface.cases
                dir = case_output_dir(interface, iteration, member, c)
                bounds = default_z_bounds(c)
                for name in interface.vars
                    var = model_scored_var(
                        dir,
                        name,
                        c,
                        pool_vars[scored_name(c, name)];
                        window = score_window(c),
                        bounds,
                    )
                    ClimaCalibrate.EnsembleBuilder.fill_g_ens_col!(builder, member, var)
                end
            end
        catch e
            @warn "Member $member failed the observation map; filling with NaN" exception =
                (e, catch_backtrace())
            ClimaCalibrate.EnsembleBuilder.fill_g_ens_col!(builder, member, NaN)
        end
    end

    ClimaCalibrate.EnsembleBuilder.is_complete(builder) || @warn(
        "G_ensemble for iteration $iteration is not fully filled."
    )
    g_ensemble = ClimaCalibrate.EnsembleBuilder.get_g_ensemble(builder)
    all(isnan, g_ensemble) && error(
        "Every member of iteration $iteration produced a NaN column.",
    )
    return g_ensemble
end

# --- SOCRATESInterface --------------------------------------------------------------------- #

"""
    SOCRATESInterface(; cases, output_dir, vars, transform, float_type, grids, run_kwargs, prune_output)

ClimaCalibrate model interface for SOCRATES EKI calibration.
"""
struct SOCRATESInterface{FT <: AbstractFloat, T, G, K} <:
       ClimaCalibrate.AbstractModelInterface
    cases::Vector{SOCRATES.SOCRATESCase}
    output_dir::String
    vars::Vector{String}
    transform::T
    float_type::Type{FT}
    grids::Vector{G}
    run_kwargs::K
    prune_output::Bool
end

function SOCRATESInterface(;
    cases = SOCRATES.all_cases(),
    output_dir::AbstractString,
    vars = collect(String, SOCRATES.SCORED_VARS),
    transform = ScoreTransform(),
    float_type::Type{<:AbstractFloat} = Float64,
    grids = nothing,
    run_kwargs = (;),
    prune_output::Bool = true,
)
    any(k -> haskey(run_kwargs, k), (:grid, :dz_min, :faces)) && error(
        "Pass the vertical grid as `grids`, not `run_kwargs`.",
    )
    cases = collect(SOCRATES.SOCRATESCase, cases)
    isempty(cases) && error("SOCRATESInterface needs at least one case")
    grids =
        isnothing(grids) ? [SOCRATES.socrates_grid(float_type, c) for c in cases] :
        collect(grids)
    length(grids) == length(cases) || error(
        "SOCRATESInterface needs one grid per case: got $(length(grids)) grids for $(length(cases)) cases.",
    )
    foreach(SOCRATES.validate, cases)
    isempty(vars) && error("SOCRATESInterface needs at least one scored variable")
    output_dir = abspath(output_dir)
    mkpath(output_dir)
    return SOCRATESInterface(
        cases,
        output_dir,
        collect(String, vars),
        transform,
        float_type,
        grids,
        run_kwargs,
        prune_output,
    )
end

case_output_dir(
    interface::SOCRATESInterface,
    iteration::Integer,
    member::Integer,
    c::SOCRATES.SOCRATESCase,
) = joinpath(
    ClimaCalibrate.path_to_ensemble_member(interface.output_dir, iteration, member),
    SOCRATES.case_name(c),
)

function case_grid(interface::SOCRATESInterface, c::SOCRATES.SOCRATESCase)
    idx = findfirst(x -> SOCRATES.case_name(x) == SOCRATES.case_name(c), interface.cases)
    isnothing(idx) && error("`$(SOCRATES.case_name(c))` is not one of this interface's cases.")
    return interface.grids[idx]
end

function run_case_for_member(
    interface::SOCRATESInterface,
    iteration::Integer,
    member::Integer,
    c::SOCRATES.SOCRATESCase,
)
    params = ClimaCalibrate.parameter_path(interface.output_dir, iteration, member)
    isfile(params) ||
        error("Sampled parameter file not found for member $member: $params")
    return SOCRATES.run_case(
        c;
        FT = interface.float_type,
        params,
        output_dir = case_output_dir(interface, iteration, member, c),
        verbose = false,
        grid = case_grid(interface, c),
        interface.run_kwargs...,
    )
end

function ClimaCalibrate.forward_model(interface::SOCRATESInterface, iteration, member)
    for c in interface.cases
        run_case_for_member(interface, iteration, member, c)
    end
    return nothing
end

ClimaCalibrate.observation_map(interface::SOCRATESInterface, iteration) =
    build_g_ensemble(interface, iteration)

ClimaCalibrate.model_interface_filepath(::SOCRATESInterface) =
    abspath(joinpath(@__DIR__, "SOCRATESCalibration.jl"))

ClimaCalibrate.experiment_dir(::SOCRATESInterface) =
    abspath(joinpath(@__DIR__, "..", ".."))

observations(interface::SOCRATESInterface) = observation_vector(
    interface.cases;
    transform = interface.transform,
    vars = interface.vars,
    grids = interface.grids,
    float_type = interface.float_type,
)

observation_series(interface::SOCRATESInterface) =
    EKP.ObservationSeries(EKP.combine_observations(observations(interface)))

# --- Prior Specification ------------------------------------------------------------------- #

const DEFAULT_PRIOR_SPEC = (
    condensation_evaporation_timescale = (1.0e2, 1e0, 1.0e8, 3.0),
    sublimation_deposition_timescale = (1.0e3, 1e1, 1.0e8, 3.0),
    rain_autoconversion_timescale = (5.0e4, 1e1, 1.0e8, 3.0),
    snow_autoconversion_timescale = (1.0e3, 1e1, 1.0e8, 3.0),
    cloud_liquid_terminal_velocity_scaling_factor = (1.0, 1 / 16, 8.0, 1.0),
    cloud_ice_terminal_velocity_scaling_factor = (1.0, 1 / 16, 8.0, 1.0),
    rain_terminal_velocity_scaling_factor = (1.0, 1 / 16, 8.0, 1.0),
    snow_terminal_velocity_scaling_factor = (1.0, 1 / 16, 8.0, 1.0),
)

function default_prior(spec = DEFAULT_PRIOR_SPEC)
    isempty(spec) && error("The prior needs at least one parameter")
    distributions = map(collect(pairs(spec))) do (name, (mean, lower, upper, σ))
        lower <= mean <= upper || error(
            "Prior mean $mean for `$name` is outside bounds ($lower, $upper).",
        )
        σ > 0 || error("Prior `unconstrained_σ` for `$name` must be positive, got $σ")
        constraint = EKP.ParameterDistributions.bounded(lower, upper)
        return EKP.ParameterDistributions.ParameterDistribution(
            EKP.ParameterDistributions.Parameterized(
                EKP.ParameterDistributions.Normal(
                    constraint.constrained_to_unconstrained(mean),
                    σ,
                ),
            ),
            constraint,
            String(name),
        )
    end
    return EKP.ParameterDistributions.combine_distributions(distributions)
end

prior_names(spec = DEFAULT_PRIOR_SPEC) = collect(String.(keys(spec)))

# --- EKP Assembly -------------------------------------------------------------------------- #

function build_ekp(
    interface::SOCRATESInterface,
    prior;
    ensemble_size::Int,
    T_stops = nothing,
    on_terminate::AbstractString = "continue",
    rng = Random.MersenneTwister(1234),
    verbose::Bool = true,
    kwargs...,
)
    ensemble_size >= 2 ||
        error("ensemble_size must be at least 2, got $ensemble_size")
    series = observation_series(interface)
    terminate_at = isnothing(T_stops) ? 1.0 : Float64(first(T_stops))
    ekp = EKP.EnsembleKalmanProcess(
        EKP.construct_initial_ensemble(rng, prior, ensemble_size),
        series,
        EKP.Inversion();
        scheduler = EKP.DataMisfitController(; terminate_at, on_terminate),
        rng,
        verbose,
        kwargs...,
    )
    check_observation_lengths(interface, ekp)
    return ekp
end

function check_observation_lengths(interface::SOCRATESInterface, ekp)
    series = EKP.get_observation_series(ekp)
    obs_len = length(EKP.get_obs(ekp))
    metadata = ClimaCalibrate.get_metadata_for_nth_iteration(series, 1)
    meta_len = sum(ClimaAnalysis.flattened_length, metadata)
    obs_len == meta_len || error(
        "Observation length ($obs_len) != metadata length ($meta_len).",
    )
    return nothing
end

function observation_block_report(interface::SOCRATESInterface)
    return map(interface.cases) do c
        ref = normalized_reference_series(
            c;
            transform = interface.transform,
            vars = interface.vars,
        )
        (;
            case = SOCRATES.case_name(c),
            length = size(ref.series, 1),
            n_times = size(ref.series, 2),
            blocks = ref.ranges,
            pool_vars = ref.pool_vars,
        )
    end
end

# --- Calibration Execution Infrastructure -------------------------------------------------- #

const DEFAULT_EMPTY_POOL_TIMEOUT = 7200

struct WorkerPoolExecutor{P <: Distributed.AbstractWorkerPool}
    pool::P
    empty_pool_timeout::Int
end

WorkerPoolExecutor(
    pool::Distributed.AbstractWorkerPool;
    empty_pool_timeout::Integer = DEFAULT_EMPTY_POOL_TIMEOUT,
) = WorkerPoolExecutor(pool, Int(empty_pool_timeout))

function run_tasks(f, tasks, executor::WorkerPoolExecutor)
    (; pool, empty_pool_timeout) = executor
    results = Vector{Any}(nothing, length(tasks))
    pending = collect(enumerate(tasks))
    inflight = Threads.Atomic{Int}(0)
    t_last_available = time()
    @sync while !isempty(pending)
        if isempty(pool.workers)
            inflight[] > 0 && (t_last_available = time())
            waited = time() - t_last_available
            inflight[] == 0 && waited > empty_pool_timeout && error(
                "No workers available for $(round(Int, waited)) s (timeout $(empty_pool_timeout) s).",
            )
            sleep(1)
            continue
        end
        t_last_available = time()
        i, task = pop!(pending)
        worker = take!(pool)
        Threads.atomic_add!(inflight, 1)
        @async try
            results[i] = Distributed.remotecall_fetch(f, worker, task)
        catch e
            @error "Task failed on worker $worker" task exception = e
        finally
            Threads.atomic_sub!(inflight, 1)
            put!(pool, worker)
        end
    end
    return results
end

# --- Calibration Driver -------------------------------------------------------------------- #

function set_field(x, name::Symbol, value)
    T = typeof(x)
    name in fieldnames(T) ||
        error("$(nameof(T)) has no field `$name`; it has $(fieldnames(T))")
    return T((f === name ? value : getfield(x, f) for f in fieldnames(T))...)
end

function apply_terminate_at(ekp, T_stops)
    isnothing(T_stops) && return ekp
    isempty(T_stops) && return ekp
    scheduler = ekp.scheduler
    hasproperty(scheduler, :terminate_at) || return ekp
    T = sum(EKP.get_Δt(ekp); init = zero(eltype(EKP.get_Δt(ekp))))
    next = findfirst(>(T), T_stops)
    isnothing(next) && return ekp
    target = T_stops[next]
    target > scheduler.terminate_at || return ekp
    @info "T_stops: advancing terminate_at" from = scheduler.terminate_at to = target T
    return set_field(ekp, :scheduler, set_field(scheduler, :terminate_at, target))
end

accumulated_T(ekp) = sum(EKP.get_Δt(ekp); init = 0.0)

function use_worker_log(dir::AbstractString)
    Distributed.myid() == 1 && return nothing
    mkpath(dir)
    path = joinpath(dir, "worker_$(Distributed.myid()).log")
    io = open(path, "w")
    Base.global_logger(Logging.SimpleLogger(io))
    @info "Logging from worker $(Distributed.myid())"
    flush(io)
    return path
end

case_marker_path(interface, iteration, member, c) =
    joinpath(case_output_dir(interface, iteration, member, c), "case_completed")

case_completed(interface, iteration, member, c) =
    isfile(case_marker_path(interface, iteration, member, c))

function mark_case_completed(interface, iteration, member, c)
    path = case_marker_path(interface, iteration, member, c)
    mkpath(dirname(path))
    write(path, "completed")
    return nothing
end

function ClimaCalibrate.Calibration.run_iteration(
    backend::ClimaCalibrate.WorkerBackend,
    interface::SOCRATESInterface,
    iteration,
    ensemble_size,
    output_dir,
)
    tasks = [
        (member, c) for member in 1:ensemble_size for c in interface.cases if
        !case_completed(interface, iteration, member, c)
    ]
    total = ensemble_size * length(interface.cases)
    @info "Iteration $iteration: $(length(tasks))/$total case-runs to do" n_workers =
        length(backend.worker_pool.workers)

    if !isempty(tasks)
        executor = WorkerPoolExecutor(
            backend.worker_pool;
            empty_pool_timeout = backend.empty_pool_timeout,
        )
        run_one = task -> begin
            member, c = task
            run_case_for_member(interface, iteration, member, c)
        end
        results = run_tasks(run_one, tasks, executor)
        for ((member, c), result) in zip(tasks, results)
            isnothing(result) || mark_case_completed(interface, iteration, member, c)
        end
    end

    failed = [
        member for member in 1:ensemble_size if
        any(c -> !case_completed(interface, iteration, member, c), interface.cases)
    ]
    for member in 1:ensemble_size
        member in failed ||
            ClimaCalibrate.write_model_completed(output_dir, iteration, member)
    end
    rate = length(failed) / ensemble_size
    rate > backend.failure_rate && error(
        "Iteration $iteration had a $(round(rate * 100; digits = 2))% failure rate.",
    )
    return nothing
end

function clear_iterations!(output_dir::AbstractString)
    isdir(output_dir) || return 0
    removed = 0
    for entry in readdir(output_dir)
        path = joinpath(output_dir, entry)
        if startswith(entry, "iteration_") && isdir(path)
            rm(path; recursive = true)
            removed += 1
        elseif entry == "eki_file.jld2"
            rm(path; force = true)
        end
    end
    return removed
end

function resume_or_initialize(ekp, prior, output_dir)
    isfile(ClimaCalibrate.ekp_path(output_dir, 1)) ||
        return ClimaCalibrate.initialize(ekp, prior, output_dir)
    next_iteration = ClimaCalibrate.last_completed_iteration(output_dir) + 1
    staged = ClimaCalibrate.load_ekp_struct(output_dir, next_iteration)
    n_staged = size(EKP.get_u_final(staged), 1)
    n_prior = EKP.ParameterDistributions.ndims(prior)
    n_staged == n_prior || error(
        "$output_dir holds a calibration over $n_staged parameters but this prior has $n_prior.",
    )
    @info "Resuming calibration" output_dir completed_iterations = next_iteration - 1
    return staged
end

"""
    calibrate(backend, ekp, interface; n_iterations, prior, T_stops, max_iter, force_termination_at_T, overwrite)

Run EKI calibration loop with staged `T_stops`.
"""
function calibrate(
    backend,
    ekp,
    interface::SOCRATESInterface;
    prior,
    n_iterations::Int,
    T_stops = nothing,
    max_iter::Int = n_iterations,
    force_termination_at_T::Bool = false,
    overwrite::Bool = false,
)
    output_dir = interface.output_dir
    if overwrite
        removed = clear_iterations!(output_dir)
        removed > 0 && @info "Discarded previous calibration" output_dir n_iterations = removed
    end
    terminal_T = isnothing(T_stops) ? Inf : Float64(last(T_stops))
    ekp = resume_or_initialize(ekp, prior, output_dir)
    ensemble_size = EKP.get_N_ens(ekp)
    @info "Running SOCRATES calibration" ensemble_size n_cases = length(interface.cases) n_iterations max_iter T_stops output_dir

    iteration = ClimaCalibrate.last_completed_iteration(output_dir) + 1
    while true
        T = accumulated_T(ekp)
        keep_going = (T < terminal_T) || (iteration <= n_iterations && !force_termination_at_T)
        keep_going || (@info "Stopping: T = $T reached final stop $terminal_T"; break)
        iteration > max_iter && (@info "Stopping: reached max_iter = $max_iter"; break)

        @info "Iteration $iteration" T
        ClimaCalibrate.Calibration.run_iteration(
            backend,
            interface,
            iteration,
            ensemble_size,
            output_dir,
        )

        ekp = ClimaCalibrate.load_ekp_struct(output_dir, iteration)
        ekp = apply_terminate_at(ekp, T_stops)
        JLD2.save_object(ClimaCalibrate.ekp_path(output_dir, iteration), ekp)
        terminate = ClimaCalibrate.observation_map_and_update!(
            ekp,
            output_dir,
            iteration,
            prior,
            interface,
        )
        ekp = ClimaCalibrate.load_ekp_struct(output_dir, iteration + 1)
        iteration += 1
        if !isnothing(terminate)
            @info "Stopping: EKP scheduler signalled termination"
            break
        end
    end
    return ekp
end

function ClimaCalibrate.analyze_iteration(
    interface::SOCRATESInterface,
    ekp,
    g_ensemble,
    prior,
    output_dir,
    iteration,
)
    @info "Iteration $iteration complete" mean_parameters =
        EKP.get_ϕ_mean_final(prior, ekp) error = last(EKP.get_error(ekp)) T = accumulated_T(ekp)
    interface.prune_output || return nothing
    removed = 0
    for member in 1:EKP.get_N_ens(ekp), c in interface.cases
        dir = case_output_dir(interface, iteration, member, c)
        isdir(dir) || continue
        for (root, _, files) in walkdir(dir), f in files
            endswith(f, ".nc") && (rm(joinpath(root, f); force = true); removed += 1)
        end
    end
    removed > 0 && @info "Pruned $removed NetCDF files from iteration $iteration"
    return nothing
end
