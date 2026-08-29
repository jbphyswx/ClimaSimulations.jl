"""
    postprocess.jl

Rerunning selected ensemble members with the full process-rate diagnostics.

A calibration writes only the scored variables, so a finished one cannot say which process drove the
result. This picks members out by their misfit and reruns them asking for
[`SOCRATES.TENDENCY_DIAGNOSTIC_VARS`](@ref).

The selection metric is the one the EKP scheduler reacts to, `Φⱼ = ½ (gⱼ - y)ᵀ Γ⁻¹ (gⱼ - y)`.
"""

"""
    member_misfits(ekp, g_ensemble)

`Φⱼ` for each column of `g_ensemble`, `NaN` for a member that failed.
"""
function member_misfits(ekp, g_ensemble::AbstractMatrix)
    y = EKP.get_obs(ekp)
    Γ = Matrix(EKP.get_obs_noise_cov(ekp))
    factorized = LinearAlgebra.cholesky(LinearAlgebra.Symmetric(Γ))
    return map(axes(g_ensemble, 2)) do j
        d = view(g_ensemble, :, j) .- y
        any(isnan, d) ? NaN : 0.5 * LinearAlgebra.dot(d, factorized \ d)
    end
end

"""Where the observation map wrote `G_ensemble` for `iteration`."""
g_ensemble_path(output_dir::AbstractString, iteration::Integer) =
    joinpath(output_dir, "iteration_" * lpad(iteration, 3, '0'), "G_ensemble.jld2")

"""
    best_members(interface; last_iteration)

The `(iteration, member, misfit)` of the lowest-misfit member over all completed iterations, and of
the lowest within the last one. Iterations with no `G_ensemble` are skipped and failed members
ignored.
"""
function best_members(
    interface::SOCRATESInterface;
    last_iteration::Integer = ClimaCalibrate.last_completed_iteration(
        interface.output_dir,
    ),
)
    last_iteration >= 1 ||
        error("No completed iterations in $(interface.output_dir); nothing to postprocess.")
    ekp = ClimaCalibrate.load_ekp_struct(interface.output_dir, last_iteration)
    candidates = NamedTuple{(:iteration, :member, :misfit), Tuple{Int, Int, Float64}}[]
    for iteration in 1:last_iteration
        path = g_ensemble_path(interface.output_dir, iteration)
        isfile(path) || continue
        for (member, φ) in enumerate(member_misfits(ekp, JLD2.load_object(path)))
            isnan(φ) || push!(candidates, (; iteration, member, misfit = φ))
        end
    end
    isempty(candidates) && error("Every ensemble member failed; there is nothing to rerun.")
    final = filter(c -> c.iteration == last_iteration, candidates)
    isempty(final) && error(
        "Every member of iteration $last_iteration failed, so there is no best final member.",
    )
    return (;
        best = argmin(c -> c.misfit, candidates),
        best_final = argmin(c -> c.misfit, final),
    )
end

"""
    rerun_member(interface, iteration, member; output_dir, period_seconds, vars, executor)

Rerun every case for one member with the full process-rate diagnostics, into
`output_dir/<case_name>`, returning the directories in `interface.cases` order.

Uses the member's own `parameters.toml` and the grid the calibration used, so this reproduces what
the calibration ran.
"""
function rerun_member(
    interface::SOCRATESInterface,
    iteration::Integer,
    member::Integer;
    output_dir::AbstractString,
    period_seconds::Real = 600,
    vars = SOCRATES.TENDENCY_DIAGNOSTIC_VARS,
    executor::SOCRATES.AbstractExecutor = SOCRATES.SerialExecutor(),
)
    params = ClimaCalibrate.parameter_path(interface.output_dir, iteration, member)
    isfile(params) ||
        error("No parameters.toml for iteration $iteration member $member: $params")
    run_one = case -> begin
        grid = case_grid(interface, case)
        SOCRATES.run_case(
            case;
            FT = interface.float_type,
            params,
            grid,
            output_dir = joinpath(output_dir, SOCRATES.case_name(case)),
            verbose = false,
            diagnostics = SOCRATES.socrates_diagnostics(
                vars;
                period_seconds,
                reduction = "average",
                n_levels = length(SOCRATES.socrates_z(grid)),
            ),
            interface.run_kwargs...,
        )
    end
    affinity_keys =
        [SOCRATES.topology_key(case_grid(interface, c)) for c in interface.cases]
    return SOCRATES.run_tasks(run_one, interface.cases, affinity_keys, executor)
end

"""
    postprocess_best_members(interface; output_dir, last_iteration, period_seconds)

Rerun the best and best-final members with full diagnostics, returning a `Dict` from label to
`(; iteration, member, misfit, dirs)`. Both are rerun even when they are the same member, so a caller
can always index either.
"""
function postprocess_best_members(
    interface::SOCRATESInterface;
    output_dir::AbstractString = joinpath(interface.output_dir, "postprocess"),
    last_iteration::Integer = ClimaCalibrate.last_completed_iteration(
        interface.output_dir,
    ),
    period_seconds::Real = 600,
    executor::SOCRATES.AbstractExecutor = SOCRATES.SerialExecutor(),
)
    selected = best_members(interface; last_iteration)
    results = Dict{String, Any}()
    for label in ("best", "best_final")
        pick = getproperty(selected, Symbol(label))
        @info "Rerunning $label member with process-rate diagnostics" pick.iteration pick.member pick.misfit
        dirs = rerun_member(
            interface,
            pick.iteration,
            pick.member;
            output_dir = joinpath(output_dir, label),
            period_seconds,
            executor,
        )
        results[label] = (; pick.iteration, pick.member, pick.misfit, dirs)
    end
    return results
end

"""
    case_budget_terms(dir, case; window, location, period)

Every term of a tendency budget written by a postprocessing run in `dir`, averaged over `window`.

Keys are budget term names: the process rates with the `location` prefix stripped, and the transport
diagnostics under their own names, together covering [`SOCRATES.MP1M_BUDGETS`](@ref) and
[`SOCRATES.TRANSPORT_BUDGETS`](@ref). Transport is registered in the grid mean only, so it is read as
such whatever `location` selects.
"""
function case_budget_terms(
    dir::AbstractString,
    case::SOCRATES.SOCRATESCase;
    window = score_window(case),
    location::AbstractString = "mp1m",
    period::AbstractString = "10m",
)
    transport = (
        "$(prefix)_$(var)" for prefix in SOCRATES.TRANSPORT_PREFIXES for
        var in SOCRATES.MP1M_BUDGET_VARS
    )
    wanted = Iterators.flatten((
        ("$(location)_$(term)" => term for term in SOCRATES.MP1M_SOURCE_TERMS),
        (name => name for name in transport),
    ))
    terms = Dict{String, Vector{Float64}}()
    z = Float64[]
    for (name, key) in wanted
        var = try
            only(values(SOCRATES.run_outputvars(dir, (name,); period)))
        catch
            continue
        end
        averaged = windowed_time_mean(var, window)
        terms[key] = collect(Float64, vec(averaged.data))
        isempty(z) &&
            (z = collect(Float64, var.dims[ClimaAnalysis.altitude_name(var)]))
    end
    isempty(terms) &&
        error("No budget terms found in $dir; was it run with `TENDENCY_DIAGNOSTIC_VARS`?")
    return terms, z
end
