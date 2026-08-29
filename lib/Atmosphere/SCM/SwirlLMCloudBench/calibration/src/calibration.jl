"""
    calibration.jl

The EKI interface for CloudBench: sampled parameters into a forward run, run output back out
as a `G` ensemble, and the prior and observations the inversion needs.

Nothing here is reachable from `SwirlLMCloudBenchSim`. A member's parameters are built by
merging the sampled file over the case's own overrides and handing the result to
`cloudbench_simulation`'s existing `params` argument, so the forward model needs no knowledge
that a calibration exists.
"""

# --- Sampled parameters into a run ---------------------------------------------------------- #

"""
    member_params(case, parameter_file; FT)

`ClimaAtmosParameters` for one ensemble member: the case's CloudBench overrides with the
sampled values layered on top.

`ClimaParams.create_toml_dict` takes either a path or a `Dict` as `override_file`, so the two
sources merge here rather than in a temporary file. The sampled values win on conflict.
"""
function member_params(
    case,
    parameter_file::AbstractString;
    FT::Type{<:AbstractFloat} = Float64,
)
    isfile(parameter_file) ||
        error("No sampled parameter file at $parameter_file.")
    experiment = Symbol(S.cloudbench_instance(case).experiment)
    overrides = SW.ClimaAtmos_SwirlLMCloudBench_toml_overrides(experiment)
    sampled = TOML.parsefile(parameter_file)
    return SW.ClimaAtmos_SwirlLMCloudBench_params(
        FT, experiment; overrides = merge(overrides, sampled),
    )
end

# --- The prior ------------------------------------------------------------------------------ #

"""
Parameter => `(mean, lower, upper, unconstrained_σ)`.

A starting point, not a prescription: pass your own `spec` to [`default_prior`](@ref). These
are the microphysics timescales a CloudBench column's condensate is most sensitive to, and the
reference runs single-moment warm rain with instantaneous condensation, which this model
cannot express — so the condensation timescale is the one parameter whose calibrated value has
a direct physical reading (how far from instantaneous the column has to sit).
"""
const DEFAULT_PRIOR_SPEC = (
    condensation_evaporation_timescale = (1.0e2, 1.0e0, 1.0e8, 3.0),
    rain_autoconversion_timescale = (5.0e4, 1.0e1, 1.0e8, 3.0),
    rain_accretion_timescale = (1.0e3, 1.0e1, 1.0e8, 3.0),
)

"""
    default_prior(spec = DEFAULT_PRIOR_SPEC)

A `ParameterDistribution` over `spec`, each entry bounded and normal in the unconstrained
space.
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

# --- Observations --------------------------------------------------------------------------- #

"""
Quantities scored by default — a starting point to edit, not the set that must be used.

`ql`, `qi`, `qr` and `qs` are height-resolved, so each costs one row per comparison level, and
each is read from a 4-D store field: `ql`/`qi` stream `q_c` and `T` together because the store
carries only total condensate.

The store publishes `lwp` alone, so `iwp`, `rwp` and `swp` are integrated from `rho` and the
condensate by their pairings' `reference_reduce`.
"""
const DEFAULT_SCORED_VARS =
    ("lwp", "olr", "asr", "cre_lw", "cre_sw", "cloud_cover", "ql", "qi", "qr", "qs", "iwp", "rwp", "swp")

"""
    scored_rows(x, name; z_src, levels)

`x` as a `(row, time)` matrix, where a scalar pairing has one row and a height-resolved one has
a row per level of `levels`, resampled from `z_src`.

Both sides of a comparison go through this, so a scored variable cannot end up reduced one way
in the reference and another in the model.

Linear in height and it refuses to extrapolate: comparison levels reaching past the source
column would otherwise be filled with invented values.
"""
function scored_rows(x, name::AbstractString; z_src, levels)
    ndims(x) == 1 && return reshape(collect(Float64, x), 1, :)
    ndims(x) == 2 || error(
        "`$name` has $(ndims(x)) dimensions; a scored variable is either one value per time \
         or height-resolved as `(z, time)`.",
    )
    isnothing(z_src) &&
        error("`$name` is height-resolved but no source height axis was given.")
    z = collect(Float64, z_src)
    size(x, 1) == length(z) || error(
        "`$name` has $(size(x, 1)) rows against $(length(z)) source levels; its leading axis \
         is not height.",
    )
    lo, hi = extrema(z)
    all(l -> lo - 1 <= l <= hi + 1, levels) || error(
        "Comparison levels $(extrema(levels)) m reach outside `$name`'s column $lo to $hi m.",
    )
    out = Matrix{Float64}(undef, length(levels), size(x, 2))
    for j in axes(x, 2)
        out[:, j] =
            CA.interp_vertical_prof(collect(Float64, levels), z, Float64.(view(x, :, j)))
    end
    return out
end

"""
    case_observation(case; vars, levels, store)

`EKP.Observation` for one case: each scored variable's time mean over the published window,
with a diagonal covariance built from the variance of that time series.

A height-resolved variable contributes one entry per level of `levels` and a scalar one entry,
so the observation length is `Σ (1 or length(levels))` rather than the number of variables.

The store carries one day at 1200 s, so the mean is over 73 samples and the variance is the
spread the column is asked to reproduce, not an instrument error.
"""
function case_observation(
    case,
    ::Type{FT} = Float64;
    vars = DEFAULT_SCORED_VARS,
    levels = nothing,
    store = CB.reference_store(case),
) where {FT <: AbstractFloat}
    isempty(vars) && error("An observation needs at least one scored variable.")
    ref = reference_comparable(case; names = collect(String, vars), store)
    z_src = ref.z
    target = isnothing(levels) ? collect(Float64, z_src) : collect(Float64, levels)
    y = FT[]
    diagonal = FT[]
    for name in vars
        rows = scored_rows(ref.data[name], name; z_src, levels = target)
        for k in axes(rows, 1)
            series = view(rows, k, :)
            push!(y, FT(Statistics.mean(series)))
            v = Statistics.var(series)
            v > 0 || error(
                "`$name`$(size(rows, 1) == 1 ? "" : " at $(target[k]) m") has zero variance \
                 in the reference, so it cannot carry a weight; drop it from `vars`.",
            )
            push!(diagonal, FT(v))
        end
    end
    return EKP.Observation(
        Dict(
            "samples" => y,
            "covariances" => LinearAlgebra.Diagonal(diagonal),
            "names" => case_name(case),
        ),
    )
end

"""
    model_altitude(run_vars)

The altitude axis shared by the height-resolved diagnostics in `run_vars`, or `nothing` when
none of them carries one.
"""
function model_altitude(run_vars, ::Type{FT} = Float64) where {FT}
    for (_, var) in run_vars
        ClimaAnalysis.has_altitude(var) || continue
        return collect(FT, var.dims[ClimaAnalysis.altitude_name(var)])
    end
    return nothing
end

"""A stable name for a case, used to label observations and member directories."""
case_name(case) = CB.cloudbench_job_id(case)

"""
    observation_vector(cases; vars)

One [`case_observation`](@ref) per case, in order.
"""
observation_vector(cases; vars = DEFAULT_SCORED_VARS, levels = nothing) =
    [case_observation(c; vars, levels) for c in cases]

# --- The interface -------------------------------------------------------------------------- #

"""
    CloudBenchInterface(; cases, output_dir, vars, float_type, grid, run_kwargs)

`ClimaCalibrate` model interface for a CloudBench EKI calibration.

`run_kwargs` are forwarded to `SwirlLMCloudBenchSim.run_case`, so `t_end`, `dt` and the model
components are the caller's to set.
"""
struct CloudBenchInterface{FT <: AbstractFloat, G, K} <:
       ClimaCalibrate.AbstractModelInterface
    cases::Vector{Any}
    output_dir::String
    vars::Vector{String}
    float_type::Type{FT}
    grid::G
    run_kwargs::K
    startup_waves::Int
    startup_wave_pause::Float64
end

function CloudBenchInterface(;
    cases,
    output_dir::AbstractString,
    vars = collect(String, DEFAULT_SCORED_VARS),
    float_type::Type{<:AbstractFloat} = Float64,
    grid = CB.cloudbench_grid(float_type; dz_min = 50.0),
    run_kwargs = (;),
    startup_waves::Integer = 0,
    startup_wave_pause::Real = 0,
)
    cases = collect(Any, cases)
    isempty(cases) && error("CloudBenchInterface needs at least one case.")
    isempty(vars) && error("CloudBenchInterface needs at least one scored variable.")
    for name in vars
        haskey(REFERENCE_PAIRINGS, name) ||
            error("`$name` has no model counterpart in `REFERENCE_PAIRINGS`.")
    end
    haskey(run_kwargs, :params) &&
        error("`params` is set per member from the sampled file; do not pass it.")
    output_dir = abspath(output_dir)
    mkpath(output_dir)
    return CloudBenchInterface(
        cases,
        output_dir,
        collect(String, vars),
        float_type,
        grid,
        run_kwargs,
        Int(startup_waves),
        Float64(startup_wave_pause),
    )
end

member_output_dir(interface::CloudBenchInterface, iteration, member, case) = joinpath(
    ClimaCalibrate.path_to_ensemble_member(interface.output_dir, iteration, member),
    case_name(case),
)

"""
    comparison_levels(interface)

The heights a height-resolved variable is compared on: the model grid's own centres, so the
model is never resampled and only the reference moves.
"""
comparison_levels(interface::CloudBenchInterface, ::Type{FT} = Float64) where {FT} =
    collect(FT, CB.cloudbench_z(interface.grid))

"""
    run_case_for_member(interface, iteration, member, case)

Run one `(member, case)` pair, returning the directory its diagnostics were written to.

This is the unit of work [`ClimaCalibrate.Calibration.run_iteration`](@ref) schedules, so it
has to be callable on a worker.
"""
function run_case_for_member(
    interface::CloudBenchInterface,
    iteration::Integer,
    member::Integer,
    case,
)
    parameter_file =
        ClimaCalibrate.parameter_path(interface.output_dir, iteration, member)
    return CB.run_case(
        case;
        FT = interface.float_type,
        output_dir = member_output_dir(interface, iteration, member, case),
        grid = interface.grid,
        params = member_params(case, parameter_file; FT = interface.float_type),
        verbose = false,
        interface.run_kwargs...,
    )
end

function ClimaCalibrate.forward_model(interface::CloudBenchInterface, iteration, member)
    for case in interface.cases
        run_case_for_member(interface, iteration, member, case)
    end
    return nothing
end

case_marker_path(interface::CloudBenchInterface, iteration, member, case) =
    joinpath(member_output_dir(interface, iteration, member, case), "case_completed")

case_completed(interface::CloudBenchInterface, iteration, member, case) =
    isfile(case_marker_path(interface, iteration, member, case))

function mark_case_completed(interface::CloudBenchInterface, iteration, member, case)
    path = case_marker_path(interface, iteration, member, case)
    mkpath(dirname(path))
    write(path, "completed")
    return nothing
end

"""
    ClimaCalibrate.Calibration.run_iteration(backend::WorkerBackend, interface::CloudBenchInterface, ...)

Run one iteration as a flat pool of `(member, case)` tasks.

With `n` sites and `m` members this is `n * m` independent tasks over the pool rather than `m`,
so the worker count is not capped by the ensemble size and a slow site cannot idle a worker
behind it. Only cases that have not already completed are scheduled, so a resumed iteration
re-runs what is missing rather than a whole member.

One grid serves every case here, so every task carries the same affinity key.
"""
function ClimaCalibrate.Calibration.run_iteration(
    backend::ClimaCalibrate.WorkerBackend,
    interface::CloudBenchInterface,
    iteration,
    ensemble_size,
    output_dir,
)
    tasks = [
        (member, case) for member in 1:ensemble_size for case in interface.cases if
        !case_completed(interface, iteration, member, case)
    ]
    total = ensemble_size * length(interface.cases)
    @info "Iteration $iteration: $(length(tasks))/$total case-runs to do" n_workers =
        length(backend.worker_pool.workers)

    if !isempty(tasks)
        executor = CB.WorkerPoolExecutor(
            backend.worker_pool;
            empty_pool_timeout = backend.empty_pool_timeout,
            interface.startup_waves,
            interface.startup_wave_pause,
        )
        affinity_keys = fill(CB.topology_key(interface.grid), length(tasks))
        run_one = task -> begin
            member, case = task
            run_case_for_member(interface, iteration, member, case)
        end
        results = CB.run_tasks(run_one, tasks, affinity_keys, executor)
        for ((member, case), result) in zip(tasks, results)
            isnothing(result) || mark_case_completed(interface, iteration, member, case)
        end
    end

    failed = [
        member for member in 1:ensemble_size if
        any(case -> !case_completed(interface, iteration, member, case), interface.cases)
    ]
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

"""
    altitude_first(var)

`var`'s data with its altitude axis leading, which is the layout [`REFERENCE_PAIRINGS`](@ref)
is written against. Data with no altitude dimension is returned unchanged.
"""
function altitude_first(var)
    data = Array(var.data)
    ClimaAnalysis.has_altitude(var) || return data
    names = collect(keys(var.dims))
    iz = findfirst(==(ClimaAnalysis.altitude_name(var)), names)
    isnothing(iz) &&
        error("Altitude is not among $(names), so its axis cannot be identified.")
    return iz == 1 ? data : permutedims(data, (iz, filter(!=(iz), 1:ndims(data))...))
end

"""
    member_g_column(interface, iteration, member)

One member's contribution to `G`: the same time means as [`case_observation`](@ref), taken
over that member's run output, stacked in case order then variable order.
"""
function member_g_column(
    interface::CloudBenchInterface,
    iteration,
    member,
    ::Type{FT} = Float64,
) where {FT <: AbstractFloat}
    return FT[value for (_, value) in member_scored_rows(interface, iteration, member, FT)]
end

"""
    member_scored_rows(interface, iteration, member) -> Vector{Pair{String, <:Real}}

One member's `G` column paired with a label per entry, in the order
[`member_g_column`](@ref) stacks them.

"""
function member_scored_rows(
    interface::CloudBenchInterface,
    iteration,
    member,
    ::Type{FT} = Float64,
) where {FT <: AbstractFloat}
    levels = comparison_levels(interface)
    needed = scored_model_vars(
        Dict(name => REFERENCE_PAIRINGS[name] for name in interface.vars),
    )
    out = Pair{String, FT}[]
    for case in interface.cases
        dir = member_output_dir(interface, iteration, member, case)
        run_vars = CB.run_outputvars(dir, needed)
        values = Dict{String, Any}(
            name => altitude_first(var) for (name, var) in run_vars
        )
        z_src = model_altitude(run_vars)
        for name in interface.vars
            rows = scored_rows(model_comparable(values, name), name; z_src, levels)
            single = size(rows, 1) == 1
            for k in axes(rows, 1)
                label =
                    single ? "$(case_name(case))/$name" :
                    "$(case_name(case))/$name@$(round(Int, levels[k]))m"
                push!(out, label => FT(Statistics.mean(view(rows, k, :))))
            end
        end
    end
    return out
end

function ClimaCalibrate.observation_map(interface::CloudBenchInterface, iteration)
    ekp = ClimaCalibrate.load_ekp_struct(interface.output_dir)
    n_members = EKP.get_N_ens(ekp)
    # From the observation itself: a height-resolved variable contributes a row per level, so
    # this is not the number of (case, variable) pairs.
    n_rows = length(EKP.get_obs(ekp))
    g_ensemble = fill(NaN, n_rows, n_members)
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

ClimaCalibrate.model_interface_filepath(::CloudBenchInterface) =
    abspath(joinpath(@__DIR__, "SwirlLMCloudBenchCalibration.jl"))

ClimaCalibrate.experiment_dir(::CloudBenchInterface) =
    abspath(joinpath(@__DIR__, ".."))

# --- EKP assembly --------------------------------------------------------------------------- #

"""
    build_ekp(interface, prior; ensemble_size, terminate_at, rng)

The `EnsembleKalmanProcess` for `interface`, with one observation per case.
"""
function build_ekp(
    interface::CloudBenchInterface,
    prior = default_prior();
    ensemble_size::Int,
    terminate_at::Real = 1.0,
    rng = Random.MersenneTwister(1234),
)
    ensemble_size > 0 || error("`ensemble_size` must be positive.")
    observations = observation_vector(interface.cases; vars = interface.vars)
    return EKP.EnsembleKalmanProcess(
        EKP.construct_initial_ensemble(rng, prior, ensemble_size),
        EKP.combine_observations(observations),
        EKP.Inversion();
        scheduler = EKP.DataMisfitController(; terminate_at),
        rng,
    )
end
