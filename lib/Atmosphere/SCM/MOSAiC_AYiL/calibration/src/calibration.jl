
#=

    The EKI calibration interface for MOSAiC AYiL

    `scoring.jl` says what is compared and how it is weighted; this turns that into
    the `(samples, covariances, names)` EKP wants, runs a member's cases, and reads
    the ensemble back.

    Every scored quantity is flattened in one order — cases outer, `vars` inner — and
    both the observation and the observation map use the same order, so the two can
    only ever line up.

=#

"""The default days a calibration runs, ascending."""
default_calibration_dates(; tops = default_calibration_tops) =
    sort!(collect(keys(tops)))

"""The default cases a calibration runs, one per day of [`default_calibration_dates`](@ref)."""
default_calibration_cases(; kwargs...) =
    [MOSAiC_AYiL.case(d) for d in default_calibration_dates(; kwargs...)]

"""
    default_MOSAiC_AYiL_grid(FT; dz_min, top)

One vertical grid for every day: the reference's own faces coarsened to `dz_min`.

The grid is identical on all 190 days, so a shared one means the ensemble compiles
and specializes once instead of once per distinct column. A day's own calibration top
still bounds what is *compared* — [`z_bounds`](@ref) — so nothing above it enters the
misfit; `top` shortens the domain itself, which only makes sense for a set of days
that share one.
"""
default_MOSAiC_AYiL_grid(
    ::Type{FT} = Float64;
    dz_min = 100,
    top = MOSAiC_AYiL.LES_TOP_FACE,
    date = first(default_calibration_dates()),
) where {FT <: AbstractFloat} = MOSAiC_AYiL.mosaic_grid(
    FT,
    MOSAiC_AYiL.case(date);
    faces = MOSAiC_AYiL.coarsen_faces_to_dz_min(
        MOSAiC_AYiL.truncate_faces_to_top(MOSAiC_AYiL.LES_FACES, top), dz_min,
    ),
)

"""
    default_calibration_grid(FT, case; dz_min)

The vertical grid a calibration runs a case on: the reference's own faces, cut at the
day's calibration top and coarsened to `dz_min`.

The top matters — above it the reference's ice is not reproducible
(`docs/design.md` §12) — and the coarsening is what makes 71 days per iteration
affordable.
"""
function default_calibration_grid(
    ::Type{FT},
    c::MOSAiC_AYiL.MOSAiCAYiLCase;
    dz_min = 50,
    tops = default_calibration_tops,
) where {FT <: AbstractFloat}
    faces = MOSAiC_AYiL.truncate_faces_to_top(
        MOSAiC_AYiL.native_faces(c), last(z_bounds(c; tops)),
    )
    return MOSAiC_AYiL.mosaic_grid(
        FT, c; faces = MOSAiC_AYiL.coarsen_faces_to_dz_min(faces, dz_min),
    )
end

"""The levels a case is scored on, from the grid it will be run on."""
grid_levels(grid, bounds) = scored_levels(MOSAiC_AYiL.mosaic_z(grid), bounds)

"""
    flat_observation(case, levels; vars, window, ...)

`(y, σ²)` for one case: every scored variable in `vars` order, in the space it is
scored in, with the diagonal of its observation covariance beside it.

`y` is not rescaled — the units live in `σ²`, which is what EKP whitens by.
"""
function flat_observation(
    c::MOSAiC_AYiL.MOSAiCAYiLCase,
    levels;
    vars = default_calibration_vars,
    window = score_window(c),
    characteristic = DEFAULT_CHARACTERISTIC,
    scaling = DEFAULT_OBS_VAR_SCALING,
    transforms = default_variable_transformations,
)
    physical = calibrated_values(c, levels, reference_fetch(c, levels); vars, window)
    y = Float64[]
    variance = Float64[]
    for name in String.(collect(vars))
        block = physical[name]
        t = get(transforms, name, NoTransform())
        v = observation_variance(
            name, block; characteristic, scaling, transform = t,
            dof = length(block), profile_dof = length(levels),
        )
        append!(y, t.(block))
        append!(variance, fill(v, length(block)))
    end
    return y, variance
end

"""
    case_observation(case, levels; vars, kwargs...)

One `EKP.Observation` per case: its scored truth and the diagonal band on it.

The covariance is diagonal by construction — the band is a statement about each
scored quantity on its own (`observation_variance`), not an estimate of how they
covary — so there is nothing here for EKP to reduce the rank of. A run that wants a
sampled, low-rank covariance instead builds it with `EKP.tsvd_cov_from_samples` and
passes it in place of this.
"""
function case_observation(
    c::MOSAiC_AYiL.MOSAiCAYiLCase,
    levels;
    vars = default_calibration_vars,
    kwargs...,
)
    y, variance = flat_observation(c, levels; vars, kwargs...)
    return EKP.Observation(
        Dict(
            "samples" => y,
            "covariances" => LinearAlgebra.Diagonal(variance),
            "names" => MOSAiC_AYiL.case_name(c),
        ),
    )
end

"""
    MOSAiCInterface(; cases, output_dir, vars, float_type, grids, run_kwargs, dz_min)

The `ClimaCalibrate` model interface for an AYiL calibration.

Holds the grid per case, because the observation is built on the levels the runs will
use: a member scored on a grid the truth was not built for would be compared against
a resampling of it.
"""
struct MOSAiCInterface{FT <: AbstractFloat, G, K} <:
       ClimaCalibrate.AbstractModelInterface
    cases::Vector{MOSAiC_AYiL.MOSAiCAYiLCase}
    output_dir::String
    vars::Vector{String}
    float_type::Type{FT}
    grids::Vector{G}
    run_kwargs::K
end

function MOSAiCInterface(;
    output_dir::AbstractString,
    dates = default_calibration_dates(),
    cases = [MOSAiC_AYiL.case(d) for d in dates],
    vars = collect(String, default_calibration_vars),
    float_type::Type{<:AbstractFloat} = Float64,
    dz_min = 100,
    grid = nothing,
    grids = isnothing(grid) ?
            [default_calibration_grid(float_type, c; dz_min) for c in cases] :
            fill(grid, length(cases)),
    run_kwargs = (;),
)
    cases = collect(MOSAiC_AYiL.MOSAiCAYiLCase, cases)
    isempty(cases) && error("MOSAiCInterface needs at least one case")
    isempty(vars) && error("MOSAiCInterface needs at least one scored variable")
    any(k -> haskey(run_kwargs, k), (:grid, :faces)) &&
        error("Pass the vertical grid as `grid` or `grids`, not `run_kwargs`.")
    length(grids) == length(cases) || error(
        "MOSAiCInterface needs one grid per case: got $(length(grids)) for \
         $(length(cases)) cases.",
    )
    output_dir = abspath(output_dir)
    mkpath(output_dir)
    return MOSAiCInterface(
        cases, output_dir, collect(String, vars), float_type, collect(grids), run_kwargs,
    )
end

function case_index(interface::MOSAiCInterface, c::MOSAiC_AYiL.MOSAiCAYiLCase)
    i = findfirst(x -> x.date == c.date, interface.cases)
    isnothing(i) && error("$(c.date) is not one of this interface's cases.")
    return i
end

case_grid(interface::MOSAiCInterface, c) = interface.grids[case_index(interface, c)]

case_levels(interface::MOSAiCInterface, c) =
    grid_levels(case_grid(interface, c), z_bounds(c))

case_output_dir(interface::MOSAiCInterface, iteration, member, c) = joinpath(
    ClimaCalibrate.path_to_ensemble_member(interface.output_dir, iteration, member),
    MOSAiC_AYiL.case_name(c),
)

"""One `EKP.Observation` per case, in the interface's case order."""
observations(interface::MOSAiCInterface) = [
    case_observation(c, case_levels(interface, c); vars = interface.vars)
    for c in interface.cases
]

"""
    observation_series(interface; minibatch_size)

The cases as an `EKP.ObservationSeries`.

Minibatched by default: 71 days per iteration is 71 forward runs per member, and EKP
draws a subset each iteration instead. A `minibatch_size` of `-1`, or the number of
cases, uses every day every iteration.
"""
function observation_series(
    interface::MOSAiCInterface;
    minibatch_size::Integer = min(8, length(interface.cases)),
)
    obs = observations(interface)
    (minibatch_size < 0 || minibatch_size >= length(obs)) &&
        return EKP.ObservationSeries(obs)
    return EKP.ObservationSeries(
        obs, EKP.RandomFixedSizeMinibatcher(minibatch_size),
    )
end

"""
    required_diagnostics(vars)

The diagnostics a run has to write for [`ClimaCalibrate.observation_map`](@ref) to
read `vars`, and nothing else.

The profiles the scored variables are, the profiles the scored paths are integrated
from, and `rhoa` for that integration. The package's own defaults are a much longer
list — radiative flux profiles, clear-sky pairs, surface fluxes — and writing those
for every member of every iteration is output nothing reads.
"""
function required_diagnostics(vars)
    profiles, paths, needed = _split_vars(vars)
    return isempty(paths) ? needed : unique(vcat(needed, "rhoa"))
end

function ClimaCalibrate.forward_model(interface::MOSAiCInterface, iteration, member)
    params = ClimaCalibrate.parameter_path(interface.output_dir, iteration, member)
    isfile(params) ||
        error("No sampled parameter file for member $member: $params")
    short_names = required_diagnostics(interface.vars)
    for c in interface.cases
        grid = case_grid(interface, c)
        MOSAiC_AYiL.run_case(
            c;
            FT = interface.float_type,
            params,
            output_dir = case_output_dir(interface, iteration, member, c),
            grid,
            diagnostics = MOSAiC_AYiL.mosaic_diagnostics(
                short_names; n_levels = length(MOSAiC_AYiL.mosaic_z(grid)),
            ),
            verbose = false,
            interface.run_kwargs...,
        )
    end
    return nothing
end

"""
    model_column(interface, iteration, member, case)

One member's scored values for one case, flattened in the same order as
[`flat_observation`](@ref).
"""
function model_column(
    interface::MOSAiCInterface,
    iteration::Integer,
    member::Integer,
    c::MOSAiC_AYiL.MOSAiCAYiLCase;
    transforms = default_variable_transformations,
)
    dir = case_output_dir(interface, iteration, member, c)
    levels = case_levels(interface, c)
    physical = calibrated_values(
        c, levels, model_fetch(dir, levels); vars = interface.vars,
        window = score_window(c),
    )
    out = Float64[]
    for name in interface.vars
        t = get(transforms, name, NoTransform())
        append!(out, t.(physical[name]))
    end
    return out
end

"""
    ClimaCalibrate.observation_map(interface, iteration)

`G_ensemble` for one iteration: a column per member over the cases EKP drew, and
`NaN` for a member whose runs did not produce one.

A `NaN` column is how EKP is told a member failed, so a member that crashed or wrote
nothing is reported rather than silently dropped.
"""
function ClimaCalibrate.observation_map(interface::MOSAiCInterface, iteration)
    ekp = ClimaCalibrate.load_ekp_struct(interface.output_dir, iteration)
    cases = minibatch_cases(interface, ekp)
    dof = sum(
        length(case_levels(interface, c)) * length(interface.vars) for c in cases
    )
    g = fill(NaN, dof, EKP.get_N_ens(ekp))
    for member in 1:EKP.get_N_ens(ekp)
        try
            g[:, member] =
                reduce(vcat, model_column(interface, iteration, member, c) for c in cases)
        catch e
            @warn "Member $member failed the observation map; its column is NaN" exception =
                (e, catch_backtrace())
        end
    end
    all(isnan, g) &&
        error("Every member of iteration $iteration produced a NaN column.")
    return g
end

"""
The cases in the minibatch EKP is comparing against this iteration.

`get_current_minibatch` gives indices into the observation vector, and
[`observations`](@ref) builds that in the interface's own case order, so the indices
are the cases.
"""
function minibatch_cases(interface::MOSAiCInterface, ekp)
    idx = EKP.get_current_minibatch(EKP.get_observation_series(ekp))
    return interface.cases[idx]
end


"""
    set_field(x, name, value)

A copy of the immutable `x` with field `name` replaced, erroring on a field it does
not have rather than silently building something else.
"""
function set_field(x, name::Symbol, value)
    T = typeof(x)
    name in fieldnames(T) ||
        error("$(nameof(T)) has no field `$name`; it has $(fieldnames(T))")
    return T((f === name ? value : getfield(x, f) for f in fieldnames(T))...)
end

"""The algorithmic time an EKP has accumulated, which is what a `T_stop` is against."""
accumulated_T(ekp) = sum(EKP.get_Δt(ekp); init = 0.0)

"""
    apply_terminate_at(ekp, T_stops)

`ekp` with its scheduler's `terminate_at` advanced to the first stop above the time
already accumulated.

The scheduler terminates the calibration when it reaches `terminate_at`, so a
schedule of stops is walked by raising that bound each time one is passed: the run
halts at each stop, which is where an ensemble can be inspected or a decision taken,
and then continues. Returns `ekp` unchanged when there is no schedule, no stop left
above the current time, or a scheduler with no such bound.
"""
function apply_terminate_at(ekp, T_stops)
    (isnothing(T_stops) || isempty(T_stops)) && return ekp
    scheduler = ekp.scheduler
    hasproperty(scheduler, :terminate_at) || return ekp
    T = accumulated_T(ekp)
    next = findfirst(>(T), T_stops)
    isnothing(next) && return ekp
    target = T_stops[next]
    target > scheduler.terminate_at || return ekp
    @info "T_stops: advancing terminate_at" from = scheduler.terminate_at to = target T
    return set_field(ekp, :scheduler, set_field(scheduler, :terminate_at, target))
end

"""
    run_calibrate(backend, ekp, interface, prior, output_dir; T_stops, max_iterations)

Run the calibration, stopping at each of `T_stops` in turn.

`ClimaCalibrate.calibrate` runs a fixed number of iterations against one termination
bound, so the loop is here instead: the bound is raised between iterations, which is
the only point at which the accumulated time is known. Ends when the last stop is
reached, when the scheduler signals termination, or at `max_iterations`.
"""
function run_calibrate(
    backend,
    ekp,
    interface::MOSAiCInterface,
    prior,
    output_dir::AbstractString;
    T_stops = nothing,
    max_iterations::Integer = 100,
)
    output_dir = abspath(output_dir)
    ensemble_size = EKP.get_N_ens(ekp)
    terminal_T = isnothing(T_stops) ? Inf : Float64(last(T_stops))
    ekp = ClimaCalibrate.initialize(ekp, prior, output_dir)
    iteration = ClimaCalibrate.last_completed_iteration(output_dir) + 1
    while true
        T = accumulated_T(ekp)
        T < terminal_T ||
            (@info "Stopping: accumulated T = $T reached the last stop $terminal_T"; break)
        iteration > max_iterations &&
            (@info "Stopping: reached max_iterations = $max_iterations"; break)

        @info "Iteration $iteration" T
        ClimaCalibrate.Calibration.run_iteration(
            backend, interface, iteration, ensemble_size, output_dir,
        )
        ekp = apply_terminate_at(
            ClimaCalibrate.load_ekp_struct(output_dir, iteration), T_stops,
        )
        JLD2.save_object(ClimaCalibrate.ekp_path(output_dir, iteration), ekp)
        terminate = ClimaCalibrate.observation_map_and_update!(
            ekp, output_dir, iteration, prior, interface,
        )
        ekp = ClimaCalibrate.load_ekp_struct(output_dir, iteration + 1)
        iteration += 1
        isnothing(terminate) ||
            (@info "Stopping: the EKP scheduler signalled termination"; break)
    end
    return ekp
end

ClimaCalibrate.model_interface_filepath(::MOSAiCInterface) =
    abspath(joinpath(@__DIR__, "MOSAiCAYiLCalibration.jl"))

ClimaCalibrate.experiment_dir(::MOSAiCInterface) =
    abspath(joinpath(@__DIR__, "..", ".."))
