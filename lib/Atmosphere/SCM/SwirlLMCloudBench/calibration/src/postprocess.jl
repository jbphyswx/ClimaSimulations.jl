"""
    postprocess.jl

Rerunning selected ensemble members with the full diagnostics.

A calibration writes only what is scored, so a finished one cannot say what drove the result.
This picks members out by their misfit and reruns them asking for
[`SwirlLMCloudBenchSim.default_diagnostic_vars`](@ref).

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

The `(iteration, member, misfit)` of the lowest-misfit member over all completed iterations, and
of the lowest within the last one. Iterations with no `G_ensemble` are skipped and failed members
ignored.
"""
function best_members(
    interface::CloudBenchInterface;
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
    rerun_member(interface, iteration, member; output_dir, vars, executor)

Rerun every case for one member with the full diagnostics, into `output_dir/<case_name>`,
returning the directories in `interface.cases` order.

Uses the member's own `parameters.toml` and the grid the calibration used, so this reproduces
what the calibration ran.
"""
function rerun_member(
    interface::CloudBenchInterface,
    iteration::Integer,
    member::Integer;
    output_dir::AbstractString,
    vars = CB.default_diagnostic_vars,
    executor::CB.AbstractExecutor = CB.SerialExecutor(),
)
    parameter_file =
        ClimaCalibrate.parameter_path(interface.output_dir, iteration, member)
    isfile(parameter_file) ||
        error("No parameters.toml for iteration $iteration member $member: $parameter_file")
    run_one = case -> CB.run_case(
        case;
        FT = interface.float_type,
        grid = interface.grid,
        params = member_params(case, parameter_file; FT = interface.float_type),
        output_dir = joinpath(output_dir, case_name(case)),
        verbose = false,
        diagnostics = CB.cloudbench_diagnostics(
            vars; n_levels = length(CB.cloudbench_z(interface.grid)),
        ),
        interface.run_kwargs...,
    )
    affinity_keys = fill(CB.topology_key(interface.grid), length(interface.cases))
    return CB.run_tasks(run_one, interface.cases, affinity_keys, executor)
end

"""
    postprocess_best_members(interface; output_dir, last_iteration, executor)

Rerun the best and best-final members with full diagnostics, returning a `Dict` from label to
`(; iteration, member, misfit, dirs)`. Both are rerun even when they are the same member, so a
caller can always index either (this is dumb).

Also this function should be refactored to specifically allow selections, best all time, best last iter, median/mean final ensemble, median/mean prior ensemble etc.
"""
function postprocess_best_members(
    interface::CloudBenchInterface;
    output_dir::AbstractString = joinpath(interface.output_dir, "postprocess"),
    last_iteration::Integer = ClimaCalibrate.last_completed_iteration(
        interface.output_dir,
    ),
    executor::CB.AbstractExecutor = CB.SerialExecutor(),
)
    selected = best_members(interface; last_iteration)
    results = Dict{String, Any}()
    for label in ("best", "best_final")
        pick = getproperty(selected, Symbol(label))
        @info "Rerunning $label member with full diagnostics" pick.iteration pick.member pick.misfit
        dirs = rerun_member(
            interface,
            pick.iteration,
            pick.member;
            output_dir = joinpath(output_dir, label),
            executor,
        )
        results[label] = (; pick.iteration, pick.member, pick.misfit, dirs)
    end
    return results
end